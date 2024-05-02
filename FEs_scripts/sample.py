import numpy as np

x = np.linspace(0,10,11)
y= np.linspace(-10,0,11)

x_grid,y_grid = np.meshgrid(x,y,indexing='ij')


print(x_grid)
print(y_grid)

mean_x = np.mean(x)
mean_y = np.mean(y)

z = np.zeros((len(x),len(y)))

# z = (x_grid-mean_x)**2/(2*(0.1**2))+ (y_grid-mean_y)**2/(2*(0.1**2))
# z = np.exp(-1*z)
z1 = np.add(x_grid,y_grid)

mask_x = (x<=5) & (x>=1)
mask_y = (y>=-5) & (y<=1)

z2 = z1[mask_x,:]
z3 = z2[:,mask_y]

print(z1)
print(z3)
