
class FDMError(Exception):
    """Exception raised for errors in the FDM."""

    def __init__(self, message="FDM error"):
        self.message = message
        super().__init__(self.message)