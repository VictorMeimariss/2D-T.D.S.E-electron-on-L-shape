import os
# AI created class to name output files if they already exist.
class FileNamer:
    def __init__(self, filepath):
        """Initialize the UniqueFileNamer with the base file path."""
        self.filepath = os.path.expanduser(filepath)
    
    def get_filename(self):
        """Generate a unique filename by appending an incrementing number if the file already exists."""
        directory, filename = os.path.split(self.filepath)
        base, ext = os.path.splitext(filename)

        # Start with _0 if the original file exists
        counter = 0
        new_filename = f"{base}_{counter}{ext}"
        new_filepath = os.path.join(directory, new_filename)

        while os.path.exists(new_filepath):
            counter += 1
            new_filename = f"{base}_{counter}{ext}"
            new_filepath = os.path.join(directory, new_filename)

        return new_filepath
