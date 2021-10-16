import matplotlib
import matplotlib.patches as patches
import numpy as np
import random

from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from numba import njit
from scipy import ndimage as ndi

matplotlib.use('Qt5Agg')

class SelectFromCollection:

    def __init__(self, ax, collection, alpha_other=0.1):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.selecting = False
        self.selected_path = None

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        print(self.fc)
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect, lineprops={"c": "b"})
        self.ind = []

    def onselect(self, verts):
        path = Path(verts,closed=True)
        print(path)
        ind = np.nonzero(path.contains_points(self.xys))[0]
        if len(ind) > 0:
            self.path = path
            self.ind = ind
            self.selecting = True
            self.fc[:, -1] = self.alpha_other
            self.fc[self.ind, -1] = 1
        # else:
        #     self.selecting = False
        #     self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

def crop_colonies(coords):

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6,6))
    fig.canvas.manager.set_window_title('Crop Colonies')
    ax = plt.subplot(1,1,1)
    ax.set_aspect('equal','box')

    pts = ax.scatter(coords[:, 0], coords[:, 1], s=3, c='r')
    selector = SelectFromCollection(ax, pts)

    def accept(event):
        if event.key == " ":
            if len(selector.ind) > 0:
                selector.selecting = False
                print("Selected points:")
                print(selector.xys[selector.ind])
                selector.fc[:, -1] = 1
                selector.collection.set_facecolors(selector.fc)
                patch = patches.PathPatch(selector.path,fc='none',lw=2,picker=True)
                ax.add_patch(patch)
                selector.fc[selector.ind, 0] = 0
                selector.fc[selector.ind, 2] = 1
                selector.collection.set_facecolors(selector.fc)
                fig.canvas.draw_idle()
        if event.key == 'escape':
            selector.selecting = False
            selector.ind = []
            selector.fc[:, -1] = 1
            selector.collection.set_facecolors(selector.fc)
            selector.canvas.draw_idle()
        if event.key == "backspace":
            if selector.selecting == False:
                print('bruh')
        if event.key == "enter":
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()

    def onpick(event):
        path = event.artist.get_path()
        print(path)
        ind = np.nonzero(path.contains_points(selector.collection.get_offsets()))[0]
        selector.fc[:, -1] = selector.alpha_other
        selector.fc[ind, -1] = 1
        selector.collection.set_facecolors(selector.fc)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", accept)
    fig.canvas.mpl_connect('pick_event', onpick)
    ax.set_title("Press enter to accept selected points.")

    plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

# fig, ax = plt.subplots()
# ax.set_title('click on points')

# line, = ax.plot(np.random.rand(100), 'o',
#                 picker=True, pickradius=5)  # 5 points tolerance

# def onpick(event):
#     thisline = event.artist
#     xdata = thisline.get_xdata()
#     ydata = thisline.get_ydata()
#     ind = event.ind
#     points = tuple(zip(xdata[ind], ydata[ind]))
#     print('onpick points:', points)

# fig.canvas.mpl_connect('pick_event', onpick)

# plt.show()