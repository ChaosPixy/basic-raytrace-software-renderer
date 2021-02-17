#include "main.h"

#include <iostream>
#include "SDL.h"
#include <vector>

struct color
{
	uint8_t r;
	uint8_t g;
	uint8_t b;

	color();

	color(const uint8_t red, const uint8_t green, const uint8_t blue)
	{
		r = red;
		g = green;
		b = blue;
	}

	bool operator== (const color& c) const
	{
		if (this->r == c.r && this->g == c.g && this->b == c.b)
			return true;

		return false;
	}

	color operator* (const double& v) const
	{
		color result(0, 0, 0);
		result.r = this->r * v;
		result.g = this->g * v;
		result.b = this->b * v;

		return result;
	}

	bool operator!= (const color& v) const
	{
		if (this->r != v.r || this->g != v.g || this->b != v.b)
			return true;

		return false;
	}

	color operator+ (const color& v) const
	{
		color result(0, 0, 0);
		result.r = this->r + v.r;
		result.g = this->g + v.g;
		result.b = this->b + v.b;

		return result;
	}
};

struct vector3
{
	double x;
	double y;
	double z;

	vector3();

	vector3(const double x, const double y, const double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	vector3 operator+ (const vector3& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x + v.x;
		result.y = this->y + v.y;
		result.z = this->z + v.z;

		return result;
	}

	vector3 operator- (const vector3& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x - v.x;
		result.y = this->y - v.y;
		result.z = this->z - v.z;

		return result;
	}

	vector3 operator* (const vector3& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x * v.x;
		result.y = this->y * v.y;
		result.z = this->z * v.z;

		return result;
	}

	vector3 operator/ (const vector3& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x / v.x;
		result.y = this->y / v.y;
		result.z = this->z / v.z;

		return result;
	}

	vector3 operator/ (const double& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x / v;
		result.y = this->y / v;
		result.z = this->z / v;

		return result;
	}

	vector3 operator* (const double& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x * v;
		result.y = this->y * v;
		result.z = this->z * v;

		return result;
	}

	vector3 operator* (const int& v) const
	{
		vector3 result(0, 0, 0);
		result.x = this->x * v;
		result.y = this->y * v;
		result.z = this->z * v;

		return result;
	}
	
	bool are_equal(const double a, const double b, const double precision) const;
	
	bool operator== (const vector3& v) const
	{
		if (are_equal(this->x, v.x, 0.00001f) && are_equal(this->y, v.y, 0.00001f) && are_equal(this->z, v.z, 0.00001f))
			return true;

		return false;
	}

	bool operator!= (const vector3& v) const
	{
		if (!are_equal(this->x, v.x, 0.00001f) || !are_equal(this->y, v.y, 0.00001f) || !are_equal(this->z, v.z, 0.00001f))
			return true;

		return false;
	}

	double length() const
	{
		return sqrt(x * x + y * y + z * z);
	}
};

struct sphere
{
	double radius;
	vector3 center;
	color col;
	double specular;
	float reflectivity;
	
	sphere();
	
	sphere(const double r, const vector3 c, const color cl, double s, float ref)
	{
		radius = r;
		center = c;
		col = cl;
		specular = s;
		reflectivity = ref;
	}

	bool operator== (const sphere& s) const
	{
		if (this->col == s.col && this->center == s.center && this->radius == s.radius)
			return true;

		return false;
	}

	bool operator!= (const sphere& s) const
	{
		if (this->col != s.col && this->center != s.center && this->radius != s.radius)
			return true;

		return false;
	}
};

enum light_type { ambient, directional, point };

struct light
{
	light();
	light(light_type t, double i, vector3 p, vector3 d);

	light_type type;
	double intensity;
	vector3 position;
	vector3 direction;
};

struct scene
{
public:
	scene();
	
	std::vector<sphere*> spheres;
	std::vector<light*> lights;
		
	~scene();
};

struct sphere_intersection
{
	sphere_intersection(sphere s, double t);
	
	sphere intersected_sphere;
	double ray_distance;
};

const color color_zero(0, 0, 0);
const vector3 vector3_zero(0, 0, 0);
const sphere sphere_zero(0, vector3_zero, color_zero, 0, 0);
const scene scene_zero = scene();

// Pointers to our window and surface
SDL_Surface* win_surface = nullptr;
SDL_Window* window = nullptr;
SDL_Renderer* renderer = nullptr; 

double canvas_width;
double canvas_height;
int max_x_coordinate;
int min_x_coordinate;
int max_y_coordinate;
int min_y_coordinate;

double view_width;
double view_height;

const double distance_to_projection_plane = 1;
const color background_color = color(0, 0, 255);

vector3 origin(0, 0, 0);

void put_pixel(const int x, const int y, const color color);
vector3 canvas_to_viewport(double x, double y);
double dot(const vector3& a, const vector3& b);
color trace_ray(const vector3 ray_source, const vector3 direction, const double ray_min, const double ray_max, uint8_t recursion_depth);
double compute_lighting(vector3 position, vector3 normal, vector3 direction_to_camera, double specular_exponent);
sphere_intersection closest_intersection(vector3 ray_source, vector3 direction, double ray_min, double ray_max);

scene scene_a = scene_zero;

void do_work();

// You must include the command line parameters for your main function to be recognized by SDL
int main(int argc, char** args)
{
	canvas_width = 600;
	canvas_height = 600;

	view_width = 1;
	view_height = 1;
	
	min_x_coordinate = canvas_width / 2 * -1;
	max_x_coordinate = canvas_width / 2;	

	min_y_coordinate = canvas_width / 2 * -1;
	max_y_coordinate = canvas_width / 2;

	
	
	// Initialize SDL. SDL_Init will return -1 if it fails.
	if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
		std::cout << "Error initializing SDL: " << SDL_GetError() << std::endl;
		system("pause");
		// End the program
		return 1;
	}

	// Create our window
	window = SDL_CreateWindow("Example", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 600, 600, SDL_WINDOW_SHOWN);

	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
	
	// Make sure creating the window succeeded
	if (!window) {
		std::cout << "Error creating window: " << SDL_GetError() << std::endl;
		system("pause");
		// End the program
		return 1;
	}

	// Get the surface from the window
	win_surface = SDL_GetWindowSurface(window);

	// Make sure getting the surface succeeded
	if (!win_surface) {
		std::cout << "Error getting surface: " << SDL_GetError() << std::endl;
		system("pause");
		// End the program
		return 1;
	}

	// Fill the window with a white rectangle
	SDL_FillRect(win_surface, nullptr, SDL_MapRGB(win_surface->format, 0, 0, 255));
	
	do_work();
		
	// Update the window display
	SDL_UpdateWindowSurface(window);
	
	// Wait
	system("pause");

	// Destroy the window. This will also destroy the surface
	SDL_DestroyWindow(window);

	// Quit SDL
	SDL_Quit();

	// End the program
	return 0;
}

void do_work()
{
	auto a = sphere(1, vector3(0, -1, 3), color(255, 0, 0), 500, 0.2f);
	scene_a.spheres.emplace_back(&a);
	auto b = sphere(1, vector3(2, 0, 4), color(0, 255, 0), 500, 0.3f);
	scene_a.spheres.emplace_back(&b);
	auto c = sphere(1, vector3(-2, 0, 4), color(255, 255, 0), 10, 0.5f);
	scene_a.spheres.emplace_back(&c);
	auto d = sphere(5000, vector3(0, -5001, 0), color(255, 0, 255), 1000, 0.4f);
	scene_a.spheres.emplace_back(&d);

	auto l1 = light(light_type::ambient, 0.2f, vector3_zero, vector3_zero);
	scene_a.lights.emplace_back(&l1);
	auto l2 = light(light_type::point, 0.6f, vector3(2, 1, 0), vector3_zero);
	scene_a.lights.emplace_back(&l2);
	auto l3 = light(light_type::directional, 0.2f, vector3_zero, vector3(1, 4, 4));
	scene_a.lights.emplace_back(&l3);

	for (auto x = min_x_coordinate; x < max_x_coordinate; x++)
	{
		for (auto y = min_y_coordinate; y < max_y_coordinate; y++)
		{
			const vector3 direction = canvas_to_viewport(x, y);
			const auto c = trace_ray(vector3_zero, direction, 1, INFINITY, 3);
			put_pixel(x, y, c);
		}
	}
		
	SDL_RenderPresent(renderer);
	
}

color::color()
{
	r = 255;
	g = 255;
	b = 255;
}

vector3::vector3()
{
	x = 0;
	y = 0;
	z = 0;
}

bool vector3::are_equal(const double a, const double b, const double precision) const
{
	return fabs(a - b) < precision;
}

sphere::sphere()
{
	radius = 0;
	center = vector3_zero;
	col = color_zero;
	
}

light::light()
{
	
}

light::light(light_type t, double i, vector3 p, vector3 d)
{
	type = t;
	intensity = i;
	position = p;
	direction = d;
}

scene::scene()
{
	spheres.reserve(3);
	lights.reserve(3);
}

scene::~scene()
{
	spheres.clear();
	lights.clear();

	for (auto* s : spheres)
		delete(s);

	for (auto* l : lights)
		delete(l);
}

sphere_intersection::sphere_intersection(sphere s, double t)
{
	intersected_sphere = s;
	ray_distance = t;
}

void put_pixel(const int x, const int y, const color color)
{
	SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
	SDL_RenderDrawPoint(renderer, canvas_width / 2 + x, canvas_height / 2 - y);
}

vector3 canvas_to_viewport(const double x, const double y)
{
	const auto result = vector3(x * view_width / canvas_width, y * view_height / canvas_height, distance_to_projection_plane);
	return result;
}

void intersect_ray_sphere(vector3 ray_source, const vector3 direction, const sphere* sphere, double &intersection1, double &intersection2)
{
	const auto r = sphere->radius;

	const auto canvas_origin = ray_source - sphere->center;

	const auto a = dot(direction, direction);
	const auto b = 2 * dot(canvas_origin, direction);
	const auto c = dot(canvas_origin, canvas_origin) - r * r;

	const auto discriminant = b * b - 4 * a * c;
	if (discriminant < 0) {
		intersection1 = INFINITY;
		intersection2 = INFINITY;
		return;
	}
	else
	{
		intersection1 = (-b + sqrt(discriminant)) / (2 * a);
		intersection2 = (-b - sqrt(discriminant)) / (2 * a);
	}
}

double dot(const vector3 &a, const vector3 &b)
{
	const double result = (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
	return result;
}

vector3 reflect_ray(vector3 normal, vector3 point)
{
	return normal * 2 * dot(normal, point) - point;
}

color trace_ray(const vector3 ray_source, const vector3 direction, const double ray_min, const double ray_max, uint8_t recursion_depth)
{
	const sphere_intersection intersection = closest_intersection(ray_source, direction, ray_min, ray_max);
	
	if (intersection.intersected_sphere == sphere_zero)
	{
		return background_color;
	}

	const vector3 position = ray_source + (direction * intersection.ray_distance);
	vector3 normal = position - intersection.intersected_sphere.center;
	normal = normal * (1.0f / normal.length());
	
	color local_color = intersection.intersected_sphere.col * compute_lighting(position, normal, direction, intersection.intersected_sphere.specular);

	float reflectivity = intersection.intersected_sphere.reflectivity;
	if (recursion_depth <= 0 || reflectivity <= 0)
	{
		return local_color;
	}

	vector3 reflected_ray = reflect_ray(normal, direction * -1);
	color reflected_color = trace_ray(position, reflected_ray, 0.001, INFINITY, recursion_depth - 1);

	return (local_color * (1 - reflectivity)) + (reflected_color * reflectivity);
}


double compute_lighting(vector3 point, vector3 normal, vector3 direction_to_camera, double specular_exponent)
{
	double i = 0.0f;
	double ray_max = INFINITY;

	for (const light* light : scene_a.lights)
	{
		if (light->type == ambient)
		{
			i += light->intensity;
		}
		else
		{
			vector3 direction = vector3_zero;
			if (light->type == light_type::point)
			{
				direction = light->position - point;
				ray_max = 1;
			}
			else
			{
				direction = light->direction;
				ray_max = INFINITY;
			}

			//Shadow			
			sphere_intersection intersection = closest_intersection(point, direction, 0.001, ray_max);			
			if (intersection.intersected_sphere != sphere_zero) continue;
			
			//Diffuse
			const double normal_dot_direction = dot(normal, direction);
			if (normal_dot_direction > 0)
			{
				i += light->intensity * normal_dot_direction / (normal.length() * point.length());
			}

			//Specular
			if (specular_exponent != -1)
			{
				vector3 reflection_vector = reflect_ray(normal, point);
				const double r_dot_v = dot(reflection_vector, direction_to_camera);

				if (r_dot_v > 0)
				{
					i += light->intensity * pow(r_dot_v / (reflection_vector.length() * direction_to_camera.length()), specular_exponent);
				}
			}
		}
	}

	if (i > 1)
		return 1;
	
	return i;
}



sphere_intersection closest_intersection(vector3 ray_source, vector3 direction, double ray_min, double ray_max)
{
	double ray_length1 = 0;
	double ray_length2 = 0;

	sphere_intersection intersection(sphere_zero, INFINITY);
	
	for (const sphere* sphere : scene_a.spheres)
	{
		intersect_ray_sphere(ray_source, direction, sphere, ray_length1, ray_length2);

		if ((ray_length1 < intersection.ray_distance) && (ray_length1 > ray_min) && (ray_length1 < ray_max))
		{
			intersection.ray_distance = ray_length1;
			intersection.intersected_sphere = *sphere;
		}
		if ((ray_length2 < intersection.ray_distance) && (ray_length2 > ray_min) && (ray_length2 < ray_max))
		{
			intersection.ray_distance = ray_length2;
			intersection.intersected_sphere = *sphere;
		}
	}	
	
	return intersection;
}