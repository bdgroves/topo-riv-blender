import requests
import os
import re


def download_google_font(font_name, output_dir="fonts", new_font_name=None):
    """Download Google font

    Parameters
    ----------
    font_name: 
    output_dir: 
    new_font_name: 
    
    Returns
    -------

    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Base URL for Google Fonts API
    url = f"https://fonts.googleapis.com/css2?family={font_name.replace(' ', '+')}"

    # Fetch the CSS file that includes the font-face rules
    response = requests.get(url)

    if response.status_code == 200:
        css_content = response.text

        # Use regex to find all URL links in the CSS content
        font_urls = re.findall(r"url\((.*?)\)", css_content)

        # Download each font file
        for font_url in font_urls:
            # Clean up the URL (remove quotes)
            font_url = font_url.strip('"').strip("'")
            if font_url.startswith("http"):
                font_response = requests.get(font_url)

                # Determine the filename
                if new_font_name:
                    # Use the new font name provided
                    font_filename = os.path.join(
                        output_dir, f"{new_font_name}.ttf"
                    )  # Change extension as needed
                else:
                    # Default to the original filename
                    font_filename = os.path.join(output_dir, font_url.split("/")[-1])

                # Save the font file
                with open(font_filename, "wb") as font_file:
                    font_file.write(font_response.content)
                print(f"Downloaded: {font_filename}")
            else:
                print(f"Invalid URL: {font_url}")  # Handle invalid URLs
    else:
        print(f"Failed to fetch font: {response.status_code}")


if __name__ == "__main__":
    label_font = snakemake.params["label_font"]
    download_google_font(label_font, new_font_name=label_font.replace(" ", "_"))
