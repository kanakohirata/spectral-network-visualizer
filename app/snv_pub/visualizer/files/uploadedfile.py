from django.conf import settings
# from django.core.files import temp as tempfile
from django.core.files.uploadedfile import UploadedFile
import os
import tempfile


class UploadedFileInTemporaryDirectory(UploadedFile):
    def __init__(self, name, content_type, size, charset, content_type_extra=None):
        _, ext = os.path.splitext(name)
        tempdir = tempfile.TemporaryDirectory(dir=settings.FILE_UPLOAD_TEMP_DIR)
        self.tempdir = tempdir
        filepath = os.path.join(tempdir.name, name)
        file = open(filepath, 'wb')
        super().__init__(file, name, content_type, size, charset, content_type_extra)

    def temporary_file_path(self):
        """Return the full path of this file."""
        return self.file.name

    def close(self):
        try:
            self.file.close()
            return self.tempdir.cleanup()
        except FileNotFoundError:
            # The file was moved or deleted before the tempfile could unlink
            # it. Still sets self.file.close_called and calls
            # self.file.file.close() before the exception.
            pass
