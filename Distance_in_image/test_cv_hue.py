
import cv2
import numpy as np

# load image with alpha channel
img = cv2.imread('IMG_0020.JPG', cv2.IMREAD_UNCHANGED)

# Convert BGR to HSV
hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
 
# define range of blue color in HSV
lower_blue = np.array([120,50,50])
upper_blue = np.array([130,255,255])
 
# Threshold the HSV image to get only blue colors
mask = cv2.inRange(hsv, lower_blue, upper_blue)
 
# Bitwise-AND mask and original image
res = cv2.bitwise_and(img, img, mask= mask)

# save output
cv2.imwrite('img_original.png', img)
cv2.imwrite('img_mask.png', mask)
cv2.imwrite('img_res.png', res)
 
cv2.imshow('frame',img)
cv2.imshow('mask',mask)
cv2.imshow('res',res)
k = cv2.waitKey(0)  


cv2.destroyAllWindows()



