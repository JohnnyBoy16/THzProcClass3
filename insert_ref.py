"""
Script to append a reference file name to a .tvl file
"""
import pdb
import os


def main():
    basedir = 'E:\\RR 2016\\THz Data\\2-4'
    filename = '2-4 Scan 12-6-2016 (res=0.25mm 100ps) - Copy.tvl'

    ref_file = 'ref 30AUG2017\\60ps waveform.txt'

    if basedir is not None:
        filename = os.path.join(basedir, filename)

    append_header(filename, ref_file)


def append_header(filename, ref_file):
    """
    """
    # at the moment this seems to delete a few very important binary characters
    # right after the header. The deletion of these characters prevents THzData
    # from functioning correctly

    # want to open the file for updating (+ symbol) and without deleting its
    # contents first see link below
    # https://docs.python.org/3/library/functions.html#open
    with open(filename, 'r+b') as fobj:
        header = list(iter(fobj.readline, b"ENDOFHEADER\r\n"))
        pdb.set_trace()
        # add the ref file name to the comment
        addition = 'COMMENT: ' + ref_file + '\r\n'
        header[1] = addition.encode('utf-8')

        # add the end of header signature
        header.append(b"ENDOFHEADER\n")

        # move back to the start of the file
        fobj.seek(0, 0)

        # write new header
        fobj.writelines(header)


if __name__ == '__main__':
    main()
