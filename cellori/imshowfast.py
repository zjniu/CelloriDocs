"""
Modification of Chris Beaumont's mpl-modest-image package to allow the use of
set_extent.
"""
from __future__ import print_function, division

import matplotlib
rcParams = matplotlib.rcParams

import matplotlib.image as mi
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
import matplotlib.cbook as cbook
import numpy as np

class ModestImage(mi.AxesImage):

    """
    Computationally modest image class.

    ModestImage is an extension of the Matplotlib AxesImage class
    better suited for the interactive display of larger images. Before
    drawing, ModestImage resamples the data array based on the screen
    resolution and view window. This has very little affect on the
    appearance of the image, but can substantially cut down on
    computation since calculations of unresolved or clipped pixels
    are skipped.

    The interface of ModestImage is the same as AxesImage. However, it
    does not currently support setting the 'extent' property. There
    may also be weird coordinate warping operations for images that
    I'm not aware of. Don't expect those to work either.
    """

    def __init__(self, stride_scale, *args, **kwargs):
        self._full_res = None
        self._sx, self._sy = None, None
        self._bounds = None
        self._scale_transform = None
        self.stride_scale = stride_scale
        super(ModestImage, self).__init__(*args, **kwargs)

    def set_data(self, A):
        """
        Set the image array

        ACCEPTS: numpy/PIL Image A
        """
        self._full_res = A
        self._A = A

        if self._A.dtype != np.uint8 and not np.can_cast(self._A.dtype,
                                                         float):
            raise TypeError("Image data can not convert to float")

        if (self._A.ndim not in (2, 3) or
                (self._A.ndim == 3 and self._A.shape[-1] not in (3, 4))):
                raise TypeError("Invalid dimensions for image data")

        self._imcache = None
        self._rgbacache = None
        self._oldxslice = None
        self._oldyslice = None
        self._sx, self._sy = None, None
        self._scale_transform = None

    def set_extent(self, extent):
        mi.AxesImage.set_extent(self, extent)
        self._scale_transform = None

    def get_array(self):
        """Override to return the full-resolution array"""
        return self._full_res

    def _get_transform(self):
        """Creates a transformation from the data limits (real extent) to the
        array limit (shape of array)."""
        if self._scale_transform is not None:
            return self._scale_transform

        if not self._extent:
            return mtransforms.IdentityTransform()

        x0 = y0 = 0.0
        y1, x1 = self._full_res.shape
        arrayLim = extent_to_bbox(x0, x1, y0, y1, self.origin)

        dataLim = self.axes.dataLim

        self._scale_transform = mtransforms.BboxTransform(dataLim, arrayLim)
        return self._scale_transform

    def _scale_to_res(self):
        """ Change self._A and _extent to render an image whose
        resolution is matched to the eventual rendering."""

        ax = self.axes
        shp = self._full_res.shape
        transform = self._get_transform()
        ss = self.stride_scale

        # Find out how we need to slice the array to make sure we match the resolution of the display.
        x0, x1, sx, y0, y1, sy = extract_matched_slices(ax, shp, transform, ss)

        # Check whether we've already calculated what we need, and if so just return without doing anything further.

        if (self._bounds is not None
            and self._sx is not None and self._sy is not None
            and sx >= self._sx and sy >= self._sy
            and x0 >= self._bounds[0] and x1 <= self._bounds[1]
            and y0 >= self._bounds[2] and y1 <= self._bounds[3]):
            return

        # Slice the array using the slices determined previously to optimally match the display
        self._A = self._full_res[y0:y1:sy, x0:x1:sx]
        self._A = cbook.safe_masked_invalid(self._A)

        # We now determine the extent of the subset of the image, by determining it first in pixel space, and converting it to the 'world' coordinates.
        extentLim = extent_to_bbox(x0, x1, y0, y1, self.origin)
        extentLim = transform.inverted().transform_bbox(extentLim)
        extent = bbox_to_extent(extentLim, self.origin)
        self.set_extent(extent)

        # Finally, we cache the current settings to avoid re-computing similar arrays in future.
        self._sx = sx
        self._sy = sy
        self._bounds = (x0, x1, y0, y1)
        self.changed()

    def draw(self, renderer, *args, **kwargs):
        self._scale_to_res()
        super(ModestImage, self).draw(renderer, *args, **kwargs)


def imshow(axes, X, stride_scale=1, cmap=None, norm=None, aspect=None,
           interpolation=None, alpha=None, vmin=None, vmax=None,
           origin=None, extent=None, shape=None, filternorm=1,
           filterrad=4.0, imlim=None, resample=None, url=None, **kwargs):
    """Similar to matplotlib's imshow command, but produces a ModestImage

    Unlike matplotlib version, must explicitly specify axes
    """

    if norm is not None:
        assert(isinstance(norm, mcolors.Normalize))
    if aspect is None:
        aspect = rcParams['image.aspect']
    axes.set_aspect(aspect)
    im = ModestImage(stride_scale, axes, cmap, norm, interpolation, origin, extent,
                     filternorm=filternorm,
                     filterrad=filterrad, resample=resample, **kwargs)

    im.set_data(X)
    im.set_alpha(alpha)
    axes._set_artist_props(im)

    if im.get_clip_path() is None:
        # image does not already have clipping set, clip to axes patch
        im.set_clip_path(axes.patch)

    # if norm is None and shape is None:
    #    im.set_clim(vmin, vmax)
    if vmin is not None or vmax is not None:
        im.set_clim(vmin, vmax)
    elif norm is None:
        im.autoscale_None()

    im.set_url(url)

    # update ax.dataLim, and, if autoscaling, set viewLim to tightly fit the image, regardless of dataLim.
    im.set_extent(im.get_extent())

    axes.images.append(im)
    im._remove_method = lambda h: axes.images.remove(h)

    return im


def extent_to_bbox(x0, x1, y0, y1, origin):
    xmin = x0
    xmax = x1
    ymin = y1 if origin == 'upper' else y0
    ymax = y0 if origin == 'upper' else y1
    corners = (xmin, ymin), (xmax, ymax)
    bbox = mtransforms.Bbox.null()
    bbox.update_from_data_xy(corners)
    return bbox


def bbox_to_extent(bbox, origin):
    x0 = bbox.xmin
    x1 = bbox.xmax
    y0 = bbox.ymax if origin == 'upper' else bbox.ymin
    y1 = bbox.ymin if origin == 'upper' else bbox.ymax
    return [x0, x1, y0, y1]


def extract_matched_slices(ax, shape, transform, ss):
    """Determine the slice parameters to use, matched to the screen.

    :param ax: Axes object to query. It's extent and pixel size
               determine the slice parameters

    :param shape: Tuple of the full image shape to slice into. Upper
               boundaries for slices will be cropped to fit within
               this shape.

    :rtype: tulpe of x0, x1, sx, y0, y1, sy

    Indexing the full resolution array as array[y0:y1:sy, x0:x1:sx] returns
    a view well-matched to the axes' resolution and extent
    """

    # Find extent in display pixels (this gives the resolution we need to sample the array to)
    ext = (ax.transAxes.transform([(1, 1)]) - ax.transAxes.transform([(0, 0)]))[0]

    # Find the extent of the axes in 'world' coordinates
    viewLim = transform.transform_bbox(ax.viewLim)
    xlim = viewLim.intervalx
    ylim = viewLim.intervaly

    dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]

    def _clip(val, hi):
        return int(max(min(val, hi), 0))

    # Determine the range of pixels to extract from the array, including a 5 pixel margin all around. We ensure that the shape of the resulting array will always be at least (1, 1) even if there is really no overlap, to avoid issues.
    y0 = _clip(min(ylim) - 5, shape[0])
    y1 = _clip(max(ylim) + 5, shape[0])
    x0 = _clip(min(xlim) - 5, shape[1])
    x1 = _clip(max(xlim) + 5, shape[1])

    # Determine the strides that can be used when extracting the array
    sy = int((max(1, ss * min((y1 - y0) / 5., np.ceil(abs(dy / ext[1]))))))
    sx = int((max(1, ss * min((x1 - x0) / 5., np.ceil(abs(dx / ext[0]))))))

    return x0, x1, sx, y0, y1, sy