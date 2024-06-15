// Globetrotter
// Minimum Cumulative Weighted Distance Computation on a Regular Grid
//
// See README.md for more information about the software.

package main

// Import modules
import (
  "bufio"
  "encoding/binary"
  "flag"
  "fmt"
  "log"
  "math"
  "os"
  "path/filepath"
  "runtime/pprof"
  "strings"
  "sync"
  "time"
)

// Set data type
type DataType = int32

// Grid describes an input grid
type Grid struct {
  x0, y0 float64
  dx, dy float64
  w, h   int
  data   []float32
  noData float32
}

// Options describes computation options
type Options struct {
  Nx, Ny int
  Move   string
  CS     string
  Roll   string
  Merge  bool
  Prefer string
}

// options describes block-specific computation options
type options struct {
  x, y   int
  w, h   int
  rx, ry bool
  wg     *sync.WaitGroup
  prefer string
}

// cell describes a computation cell
type cell struct {
  x, y int
  idx  int
  val  float32
}

// move describes a neighborhood move in number of delta cells
type move struct {
  x, y int
}

// data describes computation data
type data struct {
  sum []float32
  src []DataType
  pos []int32
  dst [][]float32
  mov []move
}

// Earth's radius in meters
const (
  EquatorialRadius = 6378137.0
  MeanRadius       = 6371008.8
  PolarRadius      = 6356752.3
  Flattening       = (EquatorialRadius - PolarRadius) / EquatorialRadius
)

// Appends cell to the heap
func push(h []cell, c cell, p []int32) []cell {
  i := int32(len(h))
  h = append(h, c)
  p[c.idx] = i
  for {
    j := (i - 1) / 2
    if h[i].val >= h[j].val {
      break
    }
    h[i], h[j] = h[j], h[i]
    p[h[i].idx], p[h[j].idx] = i, j
    i = j
  }
  return h
}

// Returns cell with the minimum value from the heap
func pop(h []cell, p []int32) ([]cell, cell) {
  i := int32(0)
  n := int32(len(h) - 1)
  c := h[0]
  h[i], h[n] = h[n], h[i]
  p[h[i].idx], p[h[n].idx] = i, -1
  val := h[i].val
  for {
    j := 2*i + 1
    if j >= n || j < 0 {
      break
    }
    if k := j + 1; k < n && h[k].val < h[j].val {
      j = k
    }
    if h[j].val >= val {
      break
    }
    h[i], h[j] = h[j], h[i]
    p[h[i].idx], p[h[j].idx] = i, j
    i = j
  }
  h = h[:n]
  return h, c
}

func update(h []cell, j int32, sum float32, pos []int32) {
  h[j].val = sum
  for {
    k := (j - 1) / 2
    if sum >= h[k].val {
      break
    }
    h[j], h[k] = h[k], h[j]
    pos[h[j].idx], pos[h[k].idx] = j, k
    j = k
  }
}

// Calculates minimum cumulative cost and source for a block of the grid
func calcCost(grid *Grid, data *data, h []cell, opts options) {
  if opts.wg != nil {
    defer opts.wg.Done()
  }

  xm, ym := opts.x+opts.w, opts.y+opts.h

  // Display block information
  log.Printf("Block [%d, %d : %d, %d], %d points\n", opts.x, opts.y, xm, ym, len(h))
  if opts.rx {
    log.Printf("Rolling enabled on the x axis\n")
  }
  if opts.ry {
    log.Printf("Rolling enabled on the y axis\n")
  }

  // Start timer
  ts := time.Now()

  // While heap is not empty
  var c cell
  for len(h) > 0 {
    // Get cell with the minimum cumulative value
    h, c = pop(h, data.pos)
    // Get cell value
    val := grid.data[c.idx]
    // For each potential move
    for i, m := range data.mov {
      // Get destination cell location
      x, y := c.x+m.x, c.y+m.y
      // Roll cell location if required
      if opts.ry {
        if y < 0 || y >= ym {
          y = -y - 1
          if y < 0 {
            y += 2 * grid.h
          }
          x = (c.x + (grid.w / 2) - m.x) % grid.w
        }
      }
      if opts.rx {
        if x < 0 || x >= xm {
          x = (x + grid.w) % grid.w
        }
      }
      // Skip cell if outside the block
      if x < opts.x || x >= xm || y < opts.y || y >= ym {
        continue
      }
      // Calculate cell index
      idx := y*grid.w + x
      // Skip if origin or no data cell
      if data.sum[idx] == 0 || data.sum[idx] == -1 {
        continue
      }
      // Calculate cumulative weighted distance
      // TODO: Check if weight constant (0.5) should be updated depending on the location
      // TODO: Proper cost calculation for >= 2nd-order neighbours (if implemented)
      sum := (grid.data[idx]+val)*data.dst[c.y][i]*0.5 + c.val
      // Skip if cell is in the heap and has less cumulative value
      j := data.pos[idx]
      if j >= 0 && h[j].val < sum {
        continue
      }
      // Skip if not minimum cumulative value
      if data.sum[idx] > 0 {
        if sum > data.sum[idx] {
          continue
        } else if sum == data.sum[idx] {
          if opts.prefer == "min" {
            if data.src[idx] < data.src[c.idx] {
              continue
            }
          } else if opts.prefer == "max" {
            if data.src[idx] > data.src[c.idx] {
              continue
            }
          } else if opts.prefer == "first" {
            continue;
          }
        }
      }
      // Set cumulative value and origin
      data.sum[idx], data.src[idx] = sum, data.src[c.idx]
      // Check if cell is in the heap
      if j >= 0 {
        // Update cell position in the heap
        update(h, j, sum, data.pos)
      } else {
        // Append cell to the heap
        h = push(h, cell{x, y, idx, sum}, data.pos)
      }
    }
  }

  // Display block information
  log.Printf("Block [%d, %d : %d, %d] finished in %s\n", opts.x, opts.y, xm, ym, time.Since(ts))
}

// Calc calculates minimum cumulative cost
func Calc(grid *Grid, src []DataType, opts Options) ([]float32, []DataType, error) {
  // Validate options
  if opts.Nx <= 0 {
    return nil, nil, fmt.Errorf("Invalid number of x blocks %d", opts.Nx)
  }
  if opts.Ny <= 0 {
    return nil, nil, fmt.Errorf("Invalid number of y blocks %d", opts.Ny)
  }
  if opts.Roll == "xy" && grid.w%2 != 0 {
    return nil, nil, fmt.Errorf("Invalid width %d (should be an even number for y-axis rolling)", grid.w)
  }

  // Create computation data
  data := &data{}

  // Set movements
  switch opts.Move {
  case "N4":
    data.mov = []move{
      {0, -1}, {-1, 0}, {1, 0}, {0, 1},
    }
  case "N8":
    data.mov = []move{
      {-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1},
    }
  default:
    return nil, nil, fmt.Errorf("Invalid movement type %s", opts.Move)
  }

  // Calculate distance lookup table
  l := len(data.mov)
  data.dst = make([][]float32, grid.h)
  switch opts.CS {
  case "planar":
    data.dst[0] = make([]float32, l)
    for i, m := range data.mov {
      x, y := float64(m.x)*grid.dx, float64(m.y)*grid.dy
      data.dst[0][i] = float32(math.Sqrt(x*x + y*y))
    }
    for j := 1; j < grid.h; j++ {
      data.dst[j] = data.dst[0]
    }
  case "spherical":
    rdy := grid.dy * math.Pi / 180
    rdx2, rdy2 := grid.dx*math.Pi/360, rdy/2
    sdx2, sdy2 := make([]float64, l), make([]float64, l)
    for i, m := range data.mov {
      sx, sy := math.Sin(float64(m.x)*rdx2), math.Sin(float64(m.y)*rdy2)
      sdx2[i], sdy2[i] = sx*sx, sy*sy
    }
    ry := (grid.y0 + grid.dy/2) * math.Pi / 180
    for j := 0; j < grid.h; j++ {
      data.dst[j] = make([]float32, l)
      cry := math.Cos(ry)
      for i, m := range data.mov {
        a := sdy2[i] + cry*math.Cos(ry+float64(m.y)*rdy)*sdx2[i]
        data.dst[j][i] = float32(MeanRadius * 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a)))
      }
      ry += rdy
    }
  case "ellipsoidal":
    rdx2 := grid.dx * math.Pi / 360
    sdx2 := make([]float64, l)
    for i, m := range data.mov {
      sx := math.Sin(float64(m.x) * rdx2)
      sdx2[i] = sx * sx
    }
    rdy := grid.dy * math.Pi / 180
    ry := (grid.y0 + grid.dy/2) * math.Pi / 180
    for j := 0; j < grid.h; j++ {
      data.dst[j] = make([]float32, l)
      b1 := math.Atan((1 - Flattening) * math.Tan(ry))
      cb1 := math.Cos(b1)
      for i, m := range data.mov {
        b2 := math.Atan((1 - Flattening) * math.Tan(ry+float64(m.y)*rdy))
        sb2 := math.Sin((b2 - b1) / 2)
        a := sb2*sb2 + cb1*math.Cos(b2)*sdx2[i]
        s := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))
        ss := math.Sin(s)
        sP2 := math.Sin((b1 + b2) / 2)
        sP2 *= sP2
        sQ2 := math.Sin((b2 - b1) / 2)
        sQ2 *= sQ2
        ss2 := math.Sin(s / 2)
        ss2 *= ss2
        X := (s - ss) * sP2 * (1 - sQ2) / (1 - ss2)
        Y := (s + ss) * (1 - sP2) * sQ2 / ss2
        data.dst[j][i] = float32(EquatorialRadius * (s - Flattening/2*(X+Y)))
      }
      ry += rdy
    }
  default:
    return nil, nil, fmt.Errorf("Invalid coordinate system %s", opts.CS)
  }

  // Create data grids
  l = len(grid.data)
  data.src = make([]DataType, l)
  data.pos = make([]int32, l)
  data.sum = make([]float32, l)
  for i := range grid.data {
    if grid.data[i] == grid.noData {
      data.sum[i] = -1
    } else {
      data.sum[i], data.pos[i] = -2, -1
    }
  }

  // Calculate block dimensions
  bw, bh := grid.w/opts.Nx, grid.h/opts.Ny
  if bw == 0 {
    return nil, nil, fmt.Errorf("Invalid number of x blocks %d", opts.Nx)
  }
  if bh == 0 {
    return nil, nil, fmt.Errorf("Invalid number of y blocks %d", opts.Ny)
  }

  // Create block heaps
  n := opts.Nx * opts.Ny
  h := make([][]cell, n)

  // Distribute sources to the heaps
  for i := range grid.data {
    if src[i] == 0 {
      continue
    }
    data.sum[i], data.src[i] = 0, src[i]
    if data.pos[i] >= 0 {
      continue
    }
    px := i % grid.w
    py := i / grid.w
    var k int
    if n == 1 {
      k = 0
    } else {
      j := px / bw
      if j == opts.Nx {
        j--
      }
      k = py / bh
      if k == opts.Ny {
        k--
      }
      k = k*opts.Nx + j
    }
    h[k] = append(h[k], cell{px, py, i, 0})
    data.pos[i] = int32(len(h[k]) - 1)
  }

  // Check if single block
  if n == 1 {
    // Calculate minimum cost
    calcCost(grid, data, h[0], options{w: grid.w, h: grid.h, rx: opts.Roll != "", ry: opts.Roll == "xy", prefer: opts.Prefer})
    // Return
    return data.sum, data.src, nil
  }

  var wg sync.WaitGroup

  // Calculate minimum cumulative values for each block concurrently
  for i, by := opts.Ny-1, grid.h%opts.Ny; i >= 0; i-- {
    y := i * bh
    for j, bx := opts.Nx-1, grid.w%opts.Nx; j >= 0; j-- {
      wg.Add(1)
      go calcCost(grid, data, h[i*opts.Nx+j], options{
        x:  j * bw,
        y:  y,
        w:  bw + bx,
        h:  bh + by,
        wg: &wg,
        prefer: opts.Prefer,
      })
      bx = 0
    }
    by = 0
  }

  // Wait for all blocks to finish
  wg.Wait()

  // Return if merge flag is not set
  if !opts.Merge {
    return data.sum, data.src, nil
  }

  for i := range h {
    h[i] = nil
  }
  g := make([]cell, 0)

  // Append visited block border cells to the heap
  for y := bh - 1; y < grid.h-1; y += bh {
    idx0 := y * grid.w
    for x := 0; x < grid.w; x++ {
      idx := idx0 + x
      if data.src[idx] > 0 {
        g = push(g, cell{x, y, idx, data.sum[idx]}, data.pos)
      }
      idx += grid.w
      if data.src[idx] > 0 {
        g = push(g, cell{x, y + 1, idx, data.sum[idx]}, data.pos)
      }
    }
  }
  for y := 0; y < grid.h; y++ {
    idx0 := y * grid.w
    for x := bw - 1; x < grid.w-1; x += bw {
      idx := idx0 + x
      if data.src[idx] > 0 && data.pos[idx] < 0 {
        g = push(g, cell{x, y, idx, data.sum[idx]}, data.pos)
      }
      idx++
      if data.src[idx] > 0 && data.pos[idx] < 0 {
        g = push(g, cell{x + 1, y, idx, data.sum[idx]}, data.pos)
      }
    }
  }
  // Append visited edge cells to the heap if x-axis rolling is enabled
  if opts.Roll != "" {
    for x := 0; x < grid.w; x += grid.w - 1 {
      for y := 0; y < grid.h; y++ {
        idx := y*grid.w + x
        if data.src[idx] > 0 {
          g = push(g, cell{x, y, idx, data.sum[idx]}, data.pos)
        }
      }
    }
  }
  // Append visited edge cells to the heap if y-axis rolling is enabled
  if opts.Roll == "xy" {
    for x := 0; x < grid.w; x++ {
      for y := 0; y < grid.h; y += grid.h - 1 {
        idx := y*grid.w + x
        if data.src[idx] > 0 && data.pos[idx] < 0 {
          g = push(g, cell{x, y, idx, data.sum[idx]}, data.pos)
        }
      }
    }
  }

  // Calculate corrected minimum cumulative values
  calcCost(grid, data, g, options{w: grid.w, h: grid.h, rx: opts.Roll != "", ry: opts.Roll == "xy", prefer: opts.Prefer})

  // Return
  return data.sum, data.src, nil
}

func saveData(name string, data interface{}, grid *Grid) error {
  dt := ""
  or := binary.LittleEndian
  f, err := os.Create(name + ".img")
  if err != nil {
    return err
  }
  defer f.Close()
  bw := bufio.NewWriter(f)
  defer bw.Flush()
  var size int;
  switch d := data.(type) {
  case []byte:
    dt = "Byte"
    _, err = bw.Write(d)
    if err != nil {
      return err
    }
    size = 1
  case []int16:
    dt = "Int16"
    bs := make([]byte, 2)
    for _, v := range d {
      or.PutUint16(bs, uint16(v))
      _, err = bw.Write(bs)
      if err != nil {
        return err
      }
    }
    size = 2
  case []int32:
    dt = "Int32"
    bs := make([]byte, 4)
    for _, v := range d {
      or.PutUint32(bs, uint32(v))
      _, err = bw.Write(bs)
      if err != nil {
        return err
      }
    }
    size = 4
  case []float32:
    dt = "Float32"
    bs := make([]byte, 4)
    for _, v := range d {
      or.PutUint32(bs, math.Float32bits(v))
      _, err = bw.Write(bs)
      if err != nil {
        return err
      }
    }
    size = 4
  default:
    return fmt.Errorf("Invalid data type %T", data)
  }
  if err = bw.Flush(); err != nil {
    return nil
  }
  if err = f.Close(); err != nil {
    return err
  }
  f, err = os.Create(name + ".vrt")
  if err != nil {
    return err
  }
  defer f.Close()
  _, err = f.WriteString(
    fmt.Sprintf("<VRTDataset rasterXSize=\"%d\" rasterYSize=\"%d\">\n", grid.w, grid.h) +
      fmt.Sprintf("<GeoTransform>%.8f, %.8f, 0.0, %.8f, 0.0, %.8f</GeoTransform>\n", grid.x0, grid.dx, grid.y0, grid.dy) +
      fmt.Sprintf("<VRTRasterBand dataType=\"%s\" band=\"1\" subClass=\"VRTRawRasterBand\">\n", dt) +
      fmt.Sprintf("<SourceFilename relativeToVRT=\"1\">%s.img</SourceFilename>\n", name) +
      "<ImageOffset>0</ImageOffset>\n" +
      fmt.Sprintf("<PixelOffset>%d</PixelOffset>\n", size) +
      fmt.Sprintf("<LineOffset>%d</LineOffset>\n", size*grid.w) +
      "<ByteOrder>LSB</ByteOrder>\n" +
      "<NoDataValue>-1.0</NoDataValue>\n" +
      "</VRTRasterBand>\n" +
      "</VRTDataset>\n")
  return err
}

func main() {
  // Define computation flags
  opts := Options{}
  flag.IntVar(&opts.Nx, "nx", 4, "number of x blocks")
  flag.IntVar(&opts.Ny, "ny", 4, "number of y blocks")
  flag.StringVar(&opts.Move, "move", "N8", "type of movement (N4, N8)")
  flag.StringVar(&opts.CS, "cs", "planar", "coordinate system (planar, spherical, ellipsoidal)")
  flag.StringVar(&opts.Roll, "roll", "", "roll axes (x, xy)")
  flag.BoolVar(&opts.Merge, "merge", true, "merge blocks")
  flag.StringVar(&opts.Prefer, "prefer", "last", "preferred id if distances are equal (min, max, first, last)")

  // Define grid flags
  grid := &Grid{}
  flag.IntVar(&grid.w, "w", 0, "grid width")
  flag.IntVar(&grid.h, "h", 0, "grid height")
  flag.Float64Var(&grid.x0, "x0", 0, "origin x coordinate (left)")
  flag.Float64Var(&grid.y0, "y0", 0, "origin y coordinate (top)")
  flag.Float64Var(&grid.dx, "dx", 1, "cell width")
  flag.Float64Var(&grid.dy, "dy", -1, "cell height")
  input := flag.String("input", "", "cost grid file name")
  source := flag.String("src", "", "source grid file name")
  out := flag.String("out", "", "output file name (set 'none' to disable output)")
  noData := flag.Float64("nodata", -1, "no-data value")

  // Define testing flags
  cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")

  // Parse flags
  flag.Parse()

  // Enable CPU profiling if required
  if *cpuprofile != "" {
    f, err := os.Create(*cpuprofile)
    if err != nil {
      log.Fatal(err)
    }
    pprof.StartCPUProfile(f)
    defer pprof.StopCPUProfile()
  }

  // Read cost grid
  if grid.w <= 0 {
    log.Fatalf("Invalid width %d", grid.w)
  }
  if grid.h <= 0 {
    log.Fatalf("Invalid height %d", grid.h)
  }
  grid.data = make([]float32, grid.w*grid.h)
  grid.noData = float32(*noData)
  fmt.Printf("Loading %s", *input)
  f, err := os.Open(*input)
  defer f.Close()
  if err != nil {
    log.Fatal(err)
  }
  if err = binary.Read(f, binary.LittleEndian, grid.data); err != nil {
    log.Fatal(err)
  }
  if err = f.Close(); err != nil {
    log.Fatal(err)
  }
  fmt.Printf("\n")

  // Read source grid
  src := make([]DataType, grid.w*grid.h)
  fmt.Printf("Loading %s", *source)
  f, err = os.Open(*source)
  defer f.Close()
  if err != nil {
    log.Fatal(err)
  }
  if err = binary.Read(f, binary.LittleEndian, src); err != nil {
    log.Fatal(err)
  }
  if err = f.Close(); err != nil {
    log.Fatal(err)
  }
  fmt.Printf("\n")

  // Start timer
  start := time.Now()

  // Calculate minimum weighted cumulative distance and origin grids
  sum, src, err := Calc(grid, src, opts)
  if err != nil {
    log.Fatal(err)
  }

  // Display time elapsed
  fmt.Printf("Time elapsed: %s\n", time.Since(start))

  // Check if output should be stored
  if *out != "none" {
    var name string
    if *out == "" {
      if *input == "" {
        name = fmt.Sprintf("%dx%d", grid.w, grid.h)
        if err = saveData(name+"_val", grid.data, grid); err != nil {
          log.Fatal(err)
        }
      } else {
        name = strings.TrimSuffix(*input, filepath.Ext(*input))
      }
      name += "_" + opts.CS + "_" + opts.Move
      if opts.Roll != "" {
        name += "_R" + opts.Roll
      }
      name = fmt.Sprintf("%s_n%dx%d", name, opts.Nx, opts.Ny)
      if !opts.Merge {
        name += "_nomerge"
      }
    } else {
      name = *out
    }
    if err = saveData(name+"_sum", sum, grid); err != nil {
      log.Fatal(err)
    }
    if err = saveData(name+"_id", src, grid); err != nil {
      log.Fatal(err)
    }
  }
}
