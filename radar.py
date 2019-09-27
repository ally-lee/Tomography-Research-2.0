import math
import numpy as np
import matplotlib.pyplot as plt


class Point(object):
	def __init__(self, x, y, root):
		self.x = x
		self.y = y
		self.root = root

	def print(self):
		return "(%f, %f)" % (self.x, self.y)

	def rotate(self, theta):
		new_x = math.cos(theta) * self.x - math.sin(theta) * self.y
		new_y = math.sin(theta) * self.x + math.cos(theta) * self.y
		return Point(new_x, new_y, self.root)


class Ellipse(object):
	def __init__(self, center_x, center_y, x_radius, y_radius):
		self.center_x = center_x
		self.center_y = center_y
		self.x_radius = x_radius
		self.y_radius = y_radius

	def print(self):
		print("Ellipse:\ncenter: (%f, %f), x-radius: %f, y-radius: %f" % (self.center_x, self.center_y, self.x_radius, self.y_radius))

	def func(self, x, y):
		if ((x - self.center_x) / self.x_radius)**2 + ((y - self.center_y) / self.y_radius)**2 <= 1:
			return 1
		return 0


class Rectangle(object):
	def __init__(self, min_x, max_x, min_y, max_y):
		self.min_x = min_x
		self.max_x = max_x
		self.min_y = min_y
		self.max_y = max_y

	def print(self):
		print("Rectangle:\n[%f, %f] x [%f, %f]" % (self.min_x, self.max_x, self.min_y, self.max_y))

	def func(self, x, y):
		if self.min_x <= x <= self.max_x and self.min_y <= y <= self.max_y:
			return 1
		return 0


class Transform(object):
	def __init__(self, theta):
		self.theta = theta

	def transform(self, shape, points, num_phi, del_phi):
		result = 0

		pt = points[0].rotate(self.theta)
		f_result = shape.func(pt.x, pt.y)
		big_result = f_result * pt.root
		result += big_result

		pt = points[num_phi - 1].rotate(self.theta)
		f_result = shape.func(pt.x, pt.y)
		big_result = f_result * pt.root
		result += big_result

		for i in range(1, num_phi - 1):
			pt = points[i].rotate(self.theta)
			f_result = shape.func(pt.x, pt.y)
			big_result = f_result * pt.root
			result += (2*big_result)

		result = del_phi / 2 * result
		return result


class Data(object):
	def __init__(self, theta, t, transform):
		self.theta = theta
		self.t = t
		self.transform = transform
		self.deriv = 0

	def print(self):
		return "%f\t%f\t%f\t%f" % (self.theta, self.t, self.transform, self.deriv)


class Reconstruction(object):
	def __init__(self, x, y, delta_t):
		self.x = x
		self.y = y
		self.delta_t = delta_t

	def approximate_transform(self, theta, data_set):
		t = math.sqrt(self.x**2 + self.y**2) + \
			math.sqrt((self.x - math.cos(theta))**2 + (self.y - math.sin(theta))**2)
		index = int((t - 1) / self.delta_t)
		#print("index: %d " % index)
		y2 = data_set[index + 1].deriv
		y1 = data_set[index].deriv
		slope = (y2 - y1) / self.delta_t
		intercept = y1 - (slope * data_set[index].t)
		approximation = (slope * t) + intercept
		return approximation

	def integral_over_theta(self, delta_theta, thetas, ts, data):
		result = self.approximate_transform(thetas[0], data[0]) + \
					self.approximate_transform(thetas[len(thetas) - 1], data[len(data) - 1])
		for i in range(1, len(thetas) - 1):
			result += (2 * self.approximate_transform(thetas[i], data[i]))
		integral = result * (delta_theta / 2)
		return integral


class Image(object):
	def __init__(self, shape, num_phi):
		self.shape = shape
		self.data = []
		self.thetas = []
		self.del_theta = 0
		self.ts = []
		self.del_t = 0
		self.phis = []
		self.num_phi = num_phi
		self.del_phi = 0
		self.points = []
		self.reconstruction = []

	def set_thetas(self, num_theta):
		self.del_theta = 2. * math.pi / float(num_theta)
		theta = 0
		self.thetas.append(theta)
		for i in range(0, num_theta):
			theta += self.del_theta
			self.thetas.append(theta)
		#print(self.thetas)
		#print(len(self.thetas))

	def set_ts(self, num_t, max_t):
		self.del_t = 2. / num_t
		t = 1.
		self.ts.append(t)
		for i in range(0, num_t):
			t += self.del_t
			self.ts.append(t)
		#print(self.ts)
		#print(len(self.ts))

	def set_phis(self):
		self.del_phi = 2. * math.pi / float(self.num_phi)
		phi = 0
		self.phis.append(phi)
		for i in range(0, self.num_phi):
			phi += self.del_phi
			self.phis.append(phi)
		#print(self.phis)
		#print(len(self.phis))
		#print(self.del_phi)

	def set_points(self):
		for i in range(len(self.ts)):
			list_of_points = []
			for j in range(len(self.phis)):
				x = 0.5 + self.ts[i] / 2 * math.cos(self.phis[j])
				y = math.sin(self.phis[j]) * math.sqrt((self.ts[i]/2)**2 - 0.25)
				root = math.sqrt((self.ts[i]/2)**2 - (math.cos(self.phis[j]))**2 * 0.25)
				list_of_points.append(Point(x, y, root))
			self.points.append(list_of_points)
	
	def print_points(self):
		for i in range(len(self.ts)):
			print(self.ts[i])
			for j in range(len(self.points[i])):
				print(self.points[i][j].print() + ", ", end="")
			print("\n")

	def acquire_data(self):
		for theta in self.thetas:
			data_pts = []
			for t in range(len(self.ts)):
				curr_ellipse = Transform(theta)
				curr_transform = curr_ellipse.transform(self.shape, self.points[t], self.num_phi, self.del_phi)
				data_pt = Data(theta, self.ts[t], curr_transform)
				data_pts.append(data_pt)
				print(data_pt.print())
			for t in range(len(data_pts) - 2):
				y1 = data_pts[t].transform
				y2 = data_pts[t+1].transform
				y3 = data_pts[t+2].transform
				new_deriv = (y3 - 2*y2 + y1) / (self.del_t**2)
				data_pts[t].deriv = new_deriv
				#print(data_pts[t].print())
			self.data.append(data_pts)

	def back_project(self, num_xy):
		y = 0.5
		for i in range(num_xy + 1):
			row = []
			x = -0.5
			for j in range(num_xy + 1):
				point = Reconstruction(x, y, self.del_t)
				row.append(point.integral_over_theta(self.del_theta, self.thetas, self.ts, self.data))
				x += (1 / num_xy)
			print(row)
			self.reconstruction.append(row)
			y -= (1 / num_xy)
		np_reconstruction = np.array(self.reconstruction)
		plt.imshow(np_reconstruction, extent=[-.5, .5, -.5, .5])
		plt.show()


#******************************DRIVER******************************#

print("This program will produce a reconstruction image of an object of your choosing. It will acquire data",
	"on the object by taking integrals over ellipses along a circle with radius 1 centered at the origin.",
	"Make sure you have matplotlib installed on your computer, and don't forget to close the matplotlib window",
	"when you're finished to terminate the program.\n")
print("To get started, I'll need some information from you:")

while True:
	try:
		shape = int(input("Ellipse (1), Circle (2), or Rectangle (3)? "))
		if shape == 1 or shape == 2 or shape == 3:
			break
		else:
			print("\nOops! You need to enter either 1, 2, or 3. Let's try this again...\n")
	except ValueError:
		print("\nOops! You need to enter either 1, 2, or 3. Let's try this again...\n")

print("\nGreat! Now I'll ask you to specify the coordinates of your object.",
	"Please keep your object within the square [-0.5, 0.5] x [-0.5, 0.5] as the image will be produced in this region.")

if shape == 1:
	while True:
		try:
			center_x = float(input("Center x-coord? "))
			center_y = float(input("Center y-coord? "))
			x_radius = float(input("X-radius? "))
			y_radius = float(input("Y-radius? "))
			if x_radius > 0 and y_radius > 0:
				break
			else:
				print("\nRemember that the radii need to be positive. Let's try this again...\n")
		except ValueError:
			print("\nOops! You need to enter a valid numbers. Let's try this again...\n")
	new_shape = Ellipse(center_x, center_y, x_radius, y_radius)
elif shape == 2:
	while True:
		try:
			center_x = float(input("Center x-coord? "))
			center_y = float(input("Center y-coord? "))
			radius = float(input("Radius? "))
			if radius > 0:
				break
			else:
				print("\nRemember that the radius needs to be positive. Let's try this again...\n")
		except ValueError:
			print("\nOops! You need to enter a valid numbers. Let's try this again...\n")
	new_shape = Ellipse(center_x, center_y, radius, radius)
elif shape == 3:
	while True:
		try:
			min_x = float(input("Min x-value? "))
			max_x = float(input("Max x-value? "))
			min_y = float(input("Min y-value? "))
			max_y = float(input("Max y-value? "))
			if min_x < max_x and min_y < max_y:
				break
			else:
				print("\nRemember that the min x-value needs to be smaller than the max x-value,",
					"and the min y-value needs to be smaller than the max y-value. Let's try this again...\n")
		except ValueError:
			print("\nOops! You need to enter a valid numbers. Let's try this again...\n")
	new_shape = Rectangle(min_x, max_x, min_y, max_y)

print("\nNow you can determine the precision by choosing the number of data points the program will produce.",
	"The recommended values will produce a clear reconstruction but will make the program take a looong time to run.",
	"Use smaller numbers for a faster runtime.")

while True:
	try:
		num_theta = int(input("Number of angles along the circle? (High res: 360 / Faster: 120) "))
		num_t = int(input("Number of ellipses for each angle? (Recommended: 600) "))
		num_phi = int(input("Number of partitions for integral approximations? (Recommended: 1000) "))
		num_xy = int(input("Number of points along each axis in reconstruction image? (High res: 1000 / Medium res: 500 / Low res: 250) "))
		if num_theta > 0 and num_t > 0 and num_phi > 0 and num_xy > 0:
			break
		else:
			print("\nRemember that these all need to be positive integers. Let's try this again...\n")
	except ValueError:
		print("\nOops! You need to enter positive integers. Let's try this again...\n")

max_t = 3

image = Image(new_shape, num_phi)
image.set_thetas(num_theta)
image.set_ts(num_t, max_t)
image.set_phis()
image.set_points()
image.acquire_data()
image.back_project(num_xy)



