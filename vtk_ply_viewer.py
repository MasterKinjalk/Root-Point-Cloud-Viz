import vtk
import os
import argparse
import sys
import numpy as np


def load_ply_file(file_path):
    """
    Load a PLY file using VTK
    
    Args:
        file_path (str): Path to the PLY file
        
    Returns:
        vtk.vtkPolyData: The loaded mesh data
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist")
    
    # Create reader for PLY files
    reader = vtk.vtkPLYReader()
    reader.SetFileName(file_path)
    reader.Update()
    
    return reader.GetOutput()


def get_mesh_information(poly_data):
    """
    Extract and return basic information about the mesh
    
    Args:
        poly_data (vtk.vtkPolyData): The mesh data
        
    Returns:
        dict: Dictionary containing mesh information
    """
    info = {
        "num_points": poly_data.GetNumberOfPoints(),
        "num_cells": poly_data.GetNumberOfCells(),
        "bounds": poly_data.GetBounds(),
        "has_points_only": poly_data.GetNumberOfCells() == 0 and poly_data.GetNumberOfPoints() > 0
    }
    
    # Check if the mesh has point data
    if poly_data.GetPointData().GetNumberOfArrays() > 0:
        info["point_arrays"] = []
        for i in range(poly_data.GetPointData().GetNumberOfArrays()):
            array = poly_data.GetPointData().GetArray(i)
            if array:
                info["point_arrays"].append({
                    "name": array.GetName(),
                    "num_components": array.GetNumberOfComponents()
                })
    
    # Check if the mesh has cell data
    if poly_data.GetCellData().GetNumberOfArrays() > 0:
        info["cell_arrays"] = []
        for i in range(poly_data.GetCellData().GetNumberOfArrays()):
            array = poly_data.GetCellData().GetArray(i)
            if array:
                info["cell_arrays"].append({
                    "name": array.GetName(),
                    "num_components": array.GetNumberOfComponents()
                })
    
    return info


def convert_point_cloud_to_glyphs(poly_data, point_size=3.0):
    """
    Convert point cloud to sphere glyphs for better visualization
    
    Args:
        poly_data (vtk.vtkPolyData): The point cloud data
        point_size (float): Size of the spheres
        
    Returns:
        vtk.vtkPolyData: The glyph polydata
    """
    # Create sphere source for glyphs
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(point_size)
    sphere.SetPhiResolution(8)
    sphere.SetThetaResolution(8)
    sphere.Update()
    
    # Create glyph3D filter
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetInputData(poly_data)
    glyph.SetScaleModeToDataScalingOff()
    glyph.SetColorModeToColorByScalar()
    glyph.Update()
    
    return glyph.GetOutput()


def process_point_cloud(poly_data, point_size=3.0, use_glyphs=True):
    """
    Process a point cloud, either by setting up appropriate rendering for points
    or by converting to glyphs
    
    Args:
        poly_data (vtk.vtkPolyData): The point cloud data
        point_size (float): Size of points or glyphs
        use_glyphs (bool): Whether to use sphere glyphs
        
    Returns:
        tuple: (vtk.vtkActor, bool) - The actor and whether vertex colors were used
    """
    # Check if we have vertex colors (RGB or RGBA)
    has_vertex_colors = False
    colors = None
    
    if poly_data.GetPointData().GetScalars():
        colors = poly_data.GetPointData().GetScalars()
        if colors.GetNumberOfComponents() >= 3:  # RGB or RGBA
            has_vertex_colors = True
    
    # If we want to use glyphs and don't have too many points
    if use_glyphs and poly_data.GetNumberOfPoints() < 100000:
        # Convert point cloud to glyphs for better visualization
        glyph_data = convert_point_cloud_to_glyphs(poly_data, point_size)
        
        # Create mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(glyph_data)
        
        if has_vertex_colors:
            mapper.SetScalarModeToUsePointFieldData()
            mapper.SelectColorArray(colors.GetName())
            mapper.SetColorModeToDirectScalars()
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        if has_vertex_colors:
            # When using glyphs, we need to ensure colors are mapped correctly
            actor.GetProperty().SetColor(1, 1, 1)  # White base color to allow the scalars to show
        else:
            # No vertex colors, use a default color
            actor.GetProperty().SetColor(0.3, 0.7, 1.0)  # Light blue
    else:
        # Standard point rendering
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly_data)
        
        if has_vertex_colors:
            mapper.SetScalarModeToUsePointData()
            mapper.SelectColorArray(colors.GetName())
            mapper.SetColorModeToDirectScalars()
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        # Set point size and enable point rendering as spheres
        actor.GetProperty().SetPointSize(point_size)
        actor.GetProperty().SetRenderPointsAsSpheres(True)
        
        if not has_vertex_colors:
            # No vertex colors, use a default color
            actor.GetProperty().SetColor(0.3, 0.7, 1.0)  # Light blue
    
    return actor, has_vertex_colors


def create_pipeline_components(poly_data, wireframe=True, edge_visibility=True, point_size=3.0):
    """
    Create the standard VTK visualization pipeline components
    
    Args:
        poly_data (vtk.vtkPolyData): The mesh data
        wireframe (bool): Whether to display the mesh as wireframe
        edge_visibility (bool): Whether to show edges
        point_size (float): Size of points for point clouds
        
    Returns:
        tuple: (mapper, actor, renderer, render_window)
    """
    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.1, 0.1, 0.1)  # Dark background for better visibility
    
    # Check if we're dealing with a point cloud (no cells)
    is_point_cloud = poly_data.GetNumberOfCells() == 0 and poly_data.GetNumberOfPoints() > 0
    
    # Handle point cloud differently
    if is_point_cloud:
        # Process the point cloud
        actor, has_colors = process_point_cloud(poly_data, point_size)
        renderer.AddActor(actor)
        
        # Create a basic mapper for return value consistency
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly_data)
    else:
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly_data)
        
        # Check for color data
        if poly_data.GetPointData().GetScalars():
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToDirectScalars()
        
        # Create actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        # Set display properties
        if wireframe:
            actor.GetProperty().SetRepresentationToWireframe()
        if edge_visibility:
            actor.GetProperty().SetEdgeVisibility(True)
            actor.GetProperty().SetLineWidth(0.5)
            
        renderer.AddActor(actor)
    
    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 1000)
    render_window.SetOffScreenRendering(1)  # Enable off-screen rendering
    
    return mapper, actor, renderer, render_window


def capture_view(render_window, file_path):
    """
    Capture the current view to an image file
    
    Args:
        render_window (vtk.vtkRenderWindow): The render window
        file_path (str): Path to save the image
    """
    # Make sure the directory exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Setup window to image filter
    window_to_image = vtk.vtkWindowToImageFilter()
    window_to_image.SetInput(render_window)
    window_to_image.ReadFrontBufferOff()
    window_to_image.Update()
    
    # Set up writer based on file extension
    _, ext = os.path.splitext(file_path)
    if ext.lower() == '.png':
        writer = vtk.vtkPNGWriter()
    elif ext.lower() == '.jpg' or ext.lower() == '.jpeg':
        writer = vtk.vtkJPEGWriter()
    else:
        writer = vtk.vtkPNGWriter()
        file_path = f"{os.path.splitext(file_path)[0]}.png"
    
    writer.SetFileName(file_path)
    writer.SetInputConnection(window_to_image.GetOutputPort())
    writer.Write()
    print(f"View saved to: {file_path}")


def setup_front_view(camera, bounds):
    """
    Set up camera for front view (looking along Z axis)
    
    Args:
        camera (vtk.vtkCamera): The camera to configure
        bounds (tuple): The mesh bounds (xmin, xmax, ymin, ymax, zmin, zmax)
    """
    x_center = (bounds[0] + bounds[1]) / 2
    y_center = (bounds[2] + bounds[3]) / 2
    z_center = (bounds[4] + bounds[5]) / 2
    
    # Calculate maximum dimension for proper scaling
    x_size = bounds[1] - bounds[0]
    y_size = bounds[3] - bounds[2]
    z_size = bounds[5] - bounds[4]
    max_dim = max(x_size, y_size, z_size)
    
    # Position camera for front view
    camera.SetPosition(x_center, y_center, z_center + max_dim * 1.5)
    camera.SetFocalPoint(x_center, y_center, z_center)
    camera.SetViewUp(0, 1, 0)


def setup_top_view(camera, bounds):
    """
    Set up camera for top view (looking down from Y axis)
    
    Args:
        camera (vtk.vtkCamera): The camera to configure
        bounds (tuple): The mesh bounds (xmin, xmax, ymin, ymax, zmin, zmax)
    """
    x_center = (bounds[0] + bounds[1]) / 2
    y_center = (bounds[2] + bounds[3]) / 2
    z_center = (bounds[4] + bounds[5]) / 2
    
    # Calculate maximum dimension for proper scaling
    x_size = bounds[1] - bounds[0]
    y_size = bounds[3] - bounds[2]
    z_size = bounds[5] - bounds[4]
    max_dim = max(x_size, y_size, z_size)
    
    # Position camera for top view
    camera.SetPosition(x_center, y_center + max_dim * 1.5, z_center)
    camera.SetFocalPoint(x_center, y_center, z_center)
    camera.SetViewUp(0, 0, -1)


def setup_side_view(camera, bounds):
    """
    Set up camera for side view (looking along X axis)
    
    Args:
        camera (vtk.vtkCamera): The camera to configure
        bounds (tuple): The mesh bounds (xmin, xmax, ymin, ymax, zmin, zmax)
    """
    x_center = (bounds[0] + bounds[1]) / 2
    y_center = (bounds[2] + bounds[3]) / 2
    z_center = (bounds[4] + bounds[5]) / 2
    
    # Calculate maximum dimension for proper scaling
    x_size = bounds[1] - bounds[0]
    y_size = bounds[3] - bounds[2]
    z_size = bounds[5] - bounds[4]
    max_dim = max(x_size, y_size, z_size)
    
    # Position camera for side view
    camera.SetPosition(x_center + max_dim * 1.5, y_center, z_center)
    camera.SetFocalPoint(x_center, y_center, z_center)
    camera.SetViewUp(0, 1, 0)


def add_coordinate_axes(renderer, bounds, scale=0.25):
    """
    Add coordinate axes to the scene
    
    Args:
        renderer (vtk.vtkRenderer): The renderer
        bounds (tuple): The mesh bounds
        scale (float): Scale factor for axes size relative to model size
    """
    # Calculate the size of the model
    x_size = bounds[1] - bounds[0]
    y_size = bounds[3] - bounds[2]
    z_size = bounds[5] - bounds[4]
    max_dim = max(x_size, y_size, z_size)
    
    # Create an axes actor
    axes = vtk.vtkAxesActor()
    
    # Scale the axes to a fraction of the model size
    axes_length = max_dim * scale
    axes.SetTotalLength(axes_length, axes_length, axes_length)
    
    # Make the axis labels larger
    axes.GetXAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    axes.GetYAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    axes.GetZAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    
    # Position the axes in the lower left corner
    axes.SetPosition(bounds[0], bounds[2], bounds[4])
    
    # Add to renderer
    renderer.AddActor(axes)


def visualize_ply(ply_path, output_dir=None, interactive=False, wireframe=True, point_size=5.0, use_glyphs=True):
    """
    Main function to visualize a PLY file and generate views
    
    Args:
        ply_path (str): Path to the PLY file
        output_dir (str, optional): Directory to save output images
        interactive (bool): Whether to show an interactive window
        wireframe (bool): Whether to display the mesh as wireframe
        point_size (float): Size of points for point clouds
        use_glyphs (bool): Whether to use sphere glyphs for point clouds
    """
    try:
        # Load the PLY file
        poly_data = load_ply_file(ply_path)
        bounds = poly_data.GetBounds()
        
        # Get mesh information
        mesh_info = get_mesh_information(poly_data)
        is_point_cloud = mesh_info["num_cells"] == 0 and mesh_info["num_points"] > 0
        
        print("\nMesh Information:")
        print(f"Number of points: {mesh_info['num_points']}")
        print(f"Number of cells: {mesh_info['num_cells']}")
        
        # Print information about point data arrays
        if "point_arrays" in mesh_info and mesh_info["point_arrays"]:
            print("\nPoint Data Arrays:")
            for array in mesh_info["point_arrays"]:
                print(f"  - {array['name']} ({array['num_components']} components)")
        
        if is_point_cloud:
            print("Data type: Point Cloud (no cells)")
            # For performance reasons, adjust glyph usage based on point count
            if mesh_info["num_points"] > 50000:
                use_glyphs = False
                print(f"Large point cloud detected ({mesh_info['num_points']} points). Disabling sphere glyphs for performance.")
            
        print(f"Bounds: [{bounds[0]:.2f}, {bounds[1]:.2f}] x [{bounds[2]:.2f}, {bounds[3]:.2f}] x [{bounds[4]:.2f}, {bounds[5]:.2f}]")
        
        # Create visualization pipeline
        mapper, actor, renderer, render_window = create_pipeline_components(
            poly_data, wireframe, edge_visibility=True, point_size=point_size
        )
        
        # Add coordinate axes
        add_coordinate_axes(renderer, bounds)
        
        # Get the camera
        camera = renderer.GetActiveCamera()
        
        # If interactive mode is enabled, create an interactor
        if interactive:
            interactor = vtk.vtkRenderWindowInteractor()
            interactor.SetRenderWindow(render_window)
            render_window.SetOffScreenRendering(0)  # Disable off-screen rendering
            
            # Add orientation widget
            axes = vtk.vtkAxesActor()
            widget = vtk.vtkOrientationMarkerWidget()
            widget.SetOrientationMarker(axes)
            widget.SetInteractor(interactor)
            widget.SetEnabled(1)
            widget.InteractiveOff()
            
            # Set up initial camera view
            setup_front_view(camera, bounds)
            renderer.ResetCamera()
            
            # Initialize and start
            interactor.Initialize()
            interactor.Start()
        else:
            # Generate and save views
            if output_dir:
                # Front view
                setup_front_view(camera, bounds)
                renderer.ResetCamera()
                render_window.Render()
                front_view_path = os.path.join(output_dir, "front_view.png")
                capture_view(render_window, front_view_path)
                
                # Top view
                setup_top_view(camera, bounds)
                renderer.ResetCamera()
                render_window.Render()
                top_view_path = os.path.join(output_dir, "top_view.png")
                capture_view(render_window, top_view_path)
                
                # Side view (additional)
                setup_side_view(camera, bounds)
                renderer.ResetCamera()
                render_window.Render()
                side_view_path = os.path.join(output_dir, "side_view.png")
                capture_view(render_window, side_view_path)
                
                print(f"\nViews successfully generated in: {output_dir}")
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


def process_batch(input_dir, output_dir, file_extension=".ply", point_size=5.0):
    """
    Process multiple PLY files in a directory
    
    Args:
        input_dir (str): Directory containing PLY files
        output_dir (str): Directory to save output images
        file_extension (str): File extension to filter for
        point_size (float): Size of points for point clouds
    """
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Input directory {input_dir} does not exist")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all PLY files in the input directory
    ply_files = [f for f in os.listdir(input_dir) if f.lower().endswith(file_extension.lower())]
    
    if not ply_files:
        print(f"No {file_extension} files found in {input_dir}")
        return
    
    print(f"Found {len(ply_files)} {file_extension} files to process")
    
    # Process each file
    for i, ply_file in enumerate(ply_files, 1):
        file_path = os.path.join(input_dir, ply_file)
        file_output_dir = os.path.join(output_dir, os.path.splitext(ply_file)[0])
        
        print(f"\n[{i}/{len(ply_files)}] Processing: {ply_file}")
        visualize_ply(file_path, file_output_dir, point_size=point_size)


def main():
    """Main entry point with argument parsing"""
    parser = argparse.ArgumentParser(description="Generate front and top views of PLY files using VTK")
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest="mode", help="Operation mode")
    
    # Single file mode
    single_parser = subparsers.add_parser("single", help="Process a single PLY file")
    single_parser.add_argument("ply_file", help="Path to the PLY file")
    single_parser.add_argument("--output_dir", "-o", default="output", help="Directory to save output images")
    single_parser.add_argument("--interactive", "-i", action="store_true", help="Show interactive visualization")
    single_parser.add_argument("--solid", "-s", action="store_true", help="Show solid mesh instead of wireframe")
    single_parser.add_argument("--point_size", "-p", type=float, default=5.0, help="Point size for point clouds")
    single_parser.add_argument("--no_glyphs", action="store_true", help="Disable sphere glyphs for point clouds")
    
    # Batch mode
    batch_parser = subparsers.add_parser("batch", help="Process multiple PLY files in a directory")
    batch_parser.add_argument("input_dir", help="Directory containing PLY files")
    batch_parser.add_argument("--output_dir", "-o", default="output", help="Directory to save output images")
    batch_parser.add_argument("--extension", "-e", default=".ply", help="File extension to filter for")
    batch_parser.add_argument("--point_size", "-p", type=float, default=5.0, help="Point size for point clouds")
    
    args = parser.parse_args()
    
    # Process based on mode
    if args.mode == "single":
        return visualize_ply(
            args.ply_file, 
            args.output_dir, 
            args.interactive, 
            not args.solid,  # wireframe is True if solid is False
            args.point_size,
            not args.no_glyphs  # use_glyphs is True if no_glyphs is False
        )
    elif args.mode == "batch":
        try:
            process_batch(args.input_dir, args.output_dir, args.extension, args.point_size)
            return 0
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            return 1
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())