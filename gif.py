

from PIL import Image

import glob



# Path to your PNGs (sorted by filename)

png_files = sorted(glob.glob("./PNGs/*.png"))
# Load images

images = [Image.open(p) for p in png_files]



# Save as GIF

images[0].save(

"outputT2.gif",

save_all=True,

append_images=images[1:],

duration=100, # milliseconds per frame

loop=0 # 0 = loop forever

)

