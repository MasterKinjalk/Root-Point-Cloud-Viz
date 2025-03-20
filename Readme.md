# VTK PLY Viewer - Documentation

## Overview

VTK PLY Viewer is a high-performance tool for generating standardized views of 3D mesh files and point clouds in PLY format. It leverages the Visualization Toolkit (VTK) to provide efficient rendering and visualization of 3D structural data.

## Features

- Generate front, top, and side views of 3D meshes and point clouds
- Process single files or batch-process entire directories
- Interactive 3D visualization mode with orientation markers
- Toggle between wireframe and solid rendering for meshes
- Sphere glyphs for better point cloud visualization
- Automatic color detection for colored point clouds
- Automated image export with configurable resolution
- Detailed mesh information reporting

## Requirements

- Python 3.10.10 (Preferred)
- VTK 9.4.1 (Latest)

## Installation

1. Make sure you have Python installed
2. Install VTK:
   ```
   pip install vtk
   ```

## Basic Usage

### Single File Processing

To process a single PLY file and generate views:

```bash
python vtk_ply_viewer.py single path/to/model.ply --output_dir path/to/output
```

This will create the following files:
- `path/to/output/front_view.png`
- `path/to/output/top_view.png`
- `path/to/output/side_view.png`

### Batch Processing

To process all PLY files in a directory:

```bash
python vtk_ply_viewer.py batch models_directory --output_dir path/to/outputs
```

This will create a separate subdirectory for each PLY file inside `path/to/outputs`.

### Interactive Mode

To visualize a 3D model interactively:

```bash
python vtk_ply_viewer.py single path/to/model.ply --interactive
```

## Command Line Options

### Global Options

- `mode`: Either `single` or `batch`

### Single Mode Options

- `ply_file`: Path to the PLY file to process
- `--output_dir`, `-o`: Directory to save output images (default: "output")
- `--interactive`, `-i`: Enable interactive visualization
- `--solid`, `-s`: Show solid mesh instead of wireframe (for mesh data)
- `--point_size`, `-p`: Set point size for point cloud visualization (default: 5.0)
- `--no_glyphs`: Disable sphere glyphs for point clouds (faster for large point clouds)

### Batch Mode Options

- `input_dir`: Directory containing PLY files to process
- `--output_dir`, `-o`: Directory to save output images (default: "output")
- `--extension`, `-e`: File extension to filter for (default: ".ply")
- `--point_size`, `-p`: Set point size for point cloud visualization (default: 5.0)

## Examples

### Generate views with solid rendering:

```bash
python vtk_ply_viewer.py single path/to/model.ply --solid
```

### Visualize a point cloud with larger points:

```bash
python vtk_ply_viewer.py single path/to/pointcloud.ply --interactive --point_size 1
```

### Faster visualization for large point clouds:

```bash
python vtk_ply_viewer.py single path/to/large_pointcloud.ply --interactive --no_glyphs
```

### Process only STL files in a directory:

```bash
python vtk_ply_viewer.py batch ./models_directory --extension .stl
```

## Troubleshooting

- If point clouds appear empty in interactive mode, try increasing the point size (e.g., `--point_size 3`)
- For large point clouds (>50,000 points), use `--no_glyphs` for better performance
- If colors aren't showing correctly, check if your PLY file includes RGB or RGBA data
- If you encounter display issues in interactive mode, make sure your system supports OpenGL rendering
- For large meshes, memory usage may be high; consider using wireframe mode
- If images appear blank, check the mesh bounds and scaling - Contact me if this issue exists