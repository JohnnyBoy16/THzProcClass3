from AScanFrame import AScanFrame
from RawCScanFrame import RawCScanFrame
from BScanFrame import BScanFrame
from InterpolatedCScanFrame import InterpolatedCScanFrame


class FrameHolder:
    """
    Class to hold the frame instances
    """
    # this class allows each frame to cause actions in another by linking them
    # through a holder class
    def __init__(self, data):

        self.a_scan_frame = AScanFrame(self, data)
        self.raw_c_scan_frame = RawCScanFrame(self, data)
        self.b_scan_frame = BScanFrame(self, data)
        self.interpolated_c_scan_frame = InterpolatedCScanFrame(self, data)
