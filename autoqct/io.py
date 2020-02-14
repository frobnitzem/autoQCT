# -*- coding: utf-8 -*-

"""IO functions for molecules."""

__all__ = [ 'read_xdr' ]

from .xdrfile import xdrfile

def read_xdr(filename):
    xdr = xdrfile(filename)
    crds = []
    for frame in xdr:
        crds.append( frame.x.copy() )

    return crds
