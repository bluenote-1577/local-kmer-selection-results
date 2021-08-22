import matplotlib.pyplot as plt
import matplotlib
plt.style.use('science')
xran = range(11,20)
plt.rcParams.update({'font.size': 14})

sync_unique_percent = [0.21, 2.1, 6.03, 15.57, 27.92, 45.76, 67.24, 81.43, 88.22]
mini_unique_percent = [6.29, 9.91,15.29, 25.03, 37, 53.99, 72.32, 83.76,89.20] 
winn_unique_percent = [6.24, 9.6,15.44, 24.92, 36.73, 53.92, 72.27, 83.75,89.20] 
sync_multiplicity = [1130, 279, 78, 23, 8.1, 3.57, 2.15, 1.66, 1.48]
mini_multiplicity = [489, 136, 42.4, 14.5, 5.9, 2.98, 1.966, 1.6, 1.459]
winn_multiplicity = [491, 138.88, 42.3, 14.57, 5.918, 2.98, 1.966, 1.6, 1.459]

plt.subplot(121)
plt.plot(xran,mini_unique_percent,'o-', label = 'Open syncmer, d = 1/5', color = 'C0')
plt.plot(xran,sync_unique_percent,'o-', label = 'Minimizer, d = 1/5', color = 'C3')
#plt.plot(xran,winn_unique_percent,'o-', label = 'Robust minimizer, d = 1/5', color = 'C4')
plt.xlabel("k")
plt.ylabel("\% unique k-mers selected")
plt.legend()


plt.subplot(122)
plt.plot(xran,mini_multiplicity,'o-', color = 'C0')
plt.plot(xran,sync_multiplicity,'o-', color = 'C3')
#plt.plot(xran,winn_multiplicity,'o-', color = 'C3')
plt.xlabel("k")
plt.ylabel("average multiplicity (number of repetitions) of selected k-mers")
plt.yscale("log")

plt.show()
