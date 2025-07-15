from osgeo import gdal
import subprocess

gdal.UseExceptions()


def resample(upsample, method, input_bands, output_bands):
    for i in range(0, len(input_bands)):
        # Open the image
        dataset = gdal.Open(input_bands[i])

        # Get the raster size
        width = dataset.RasterXSize
        height = dataset.RasterYSize

        # Calculate new dimensions
        new_width = width * upsample
        new_height = height * upsample

        # Print new dimensions
        print(f"New dimensions: {new_width} x {new_height}")

        # Close the dataset
        dataset = None

        # Now you can call gdalwarp with new_width and new_height
        # Example command (this won't execute; it's just for demonstration)
        # gdalwarp -r cubic -ts {new_width} {new_height} input.tif output.tif

        # Construct the gdalwarp command
        gdalwarp_command = [
            "gdalwarp",  # Command
            "-r",
            method,  # Resampling method
            "-ts",
            str(new_width),
            str(new_height),  # New dimensions
            input_bands[i],  # Input file
            output_bands[i],  # Output file
        ]

        # Run the gdalwarp command
        try:
            subprocess.run(gdalwarp_command, check=True)
            print(
                f"Successfully created {output_bands[i]} with dimensions {new_width} x {new_height}."
            )
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running gdalwarp: {e}")


if __name__ == "__main__":
    upsample = snakemake.params["upsample"]
    method = snakemake.params["method"]
    input_bands = snakemake.input["input_bands"]
    output_bands = snakemake.output["output_bands"]

    resample(upsample, method, input_bands, output_bands)
