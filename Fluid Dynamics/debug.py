import numpy as np
import matplotlib.pyplot as plt
import imageio

w, h = 126, 126;

solid = np.ones ((h, w));
# fluid = np.ones ((h, w));

u = np.zeros ((h, w + 1));
v = np.zeros ((h + 1, w));

solid [[0, -1], :] = solid [:, [0, -1]] = -1;
# fluid [25:-1, 1:-1] = -1;
#
# u [30:40, 30] = -50;
# u [30:40, 70] = +50;

fluid = imageio.imread ('heart.png')[..., 0].astype(np.float);
# fluid = imageio.imread ('output\\108.png').astype(np.float);
fluid = ((fluid / 255) - 0.5) * 2;

print (fluid.shape);

np.savetxt ('solid.txt', solid, header = '{} {} 0.5 0.5'.format(w, h));
np.savetxt ('fluid.txt', fluid, header = '{} {} 0.5 0.5'.format(w, h));

np.savetxt ('u_test.txt', u, header = '{} {} 0. 0.5'.format (w + 1, h));
np.savetxt ('v_test.txt', v, header = '{} {} 0.5 0.'.format (w, h + 1));

u = np.loadtxt ('debug.txt');


fig, ax = plt.subplots ();

# ax.set_xticks(np.round(np.linspace(0.5, w - 0.5, num = w - 1), 1))
# ax.set_yticks(np.round(np.linspace(0.5, h - 0.5, num = h - 1), 1))
# ax.set_xticklabels([]);
# ax.set_yticklabels([]);

ax.grid();

img = ax.imshow (u);
plt.colorbar(img);
plt.show ();
