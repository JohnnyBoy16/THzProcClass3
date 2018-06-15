from THzProc.AScanFrame import AScanFrame
from THzProc.RawCScanFrame import RawCScanFrame
from THzProc.BScanFrame import BScanFrame
from THzProc.InterpolatedCScanFrame import InterpolatedCScanFrame


class FrameHolder:
    """
    Class to hold the frame instances
    """
    # this class allows each frame to cause actions in another by linking them
    # through a holder class

    def __init__(self, data, a_scan_only=True):
        """
        Constructor method. Sets all of the instance attributes and allows
        the frames to be linked together
        :param data: An instance of the THzData class
        :param a_scan_only: Either True or False depending on whether or not
            the user wants to see the frequency domain information
            (Default: False)
        """

        self.a_scan_frame = AScanFrame(self, data, a_scan_only)
        self.raw_c_scan_frame = RawCScanFrame(self, data)
        self.b_scan_frame = BScanFrame(self, data)
        self.interpolated_c_scan_frame = InterpolatedCScanFrame(self, data)
