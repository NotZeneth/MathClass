#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"
#include <iostream>

#include <algorithm>

// Définition d’un mur (segment 2D)
struct Segment {
    glm::vec2 a, b;
    Segment(glm::vec2 a_, glm::vec2 b_) : a(a_), b(b_) {}
};

struct Circle {
    glm::vec2 center;
    float radius;

    Circle(glm::vec2 c, float r) : center(c), radius(r) {}
};


bool intersect_segments(glm::vec2 p1, glm::vec2 p2, glm::vec2 q1, glm::vec2 q2, glm::vec2& intersection, float& t_out)
{
    glm::vec2 r = p2 - p1;
    glm::vec2 s = q2 - q1;
    glm::vec2 b = q1 - p1;

    glm::mat2 A(r, -s); // colonne 0 = r, colonne 1 = -s

    float det = glm::determinant(A);
    if (det == 0.f) {
        return false; // si paral
    }

    glm::vec2 x = glm::inverse(A) * b; // Forme x = [t, u]
    float t = x.x;
    float u = x.y;

    if (t >= 0.f && t <= 1.f && u >= 0.f && u <= 1.f) {
        intersection = p1 + t * r;
        t_out = t;
        return true;
    }

    return false;
}


bool intersect_segment_circle(glm::vec2 p0, glm::vec2 p1, glm::vec2 center, float radius, glm::vec2& intersection, float& t_out)
{
    glm::vec2 d = p1 - p0; // direction du déplacement
    glm::vec2 o = p0; // point de départ
    glm::vec2 c = center;

    glm::vec2 oc = o - c;

    float a = glm::dot(d, d);
    float b = 2.f * glm::dot(oc, d);
    float c_val = glm::dot(oc, oc) - radius * radius;

    float discriminant = b * b - 4.f * a * c_val;
    if (discriminant < 0.f)
        return false; // pas d'intersection

    discriminant = std::sqrt(discriminant);

    float t1 = (-b - discriminant) / (2.f * a);
    float t2 = (-b + discriminant) / (2.f * a);

    // Check si une solution est sur segment [0,1]
    if (t1 >= 0.f && t1 <= 1.f) {
        t_out = t1;
        intersection = o + t1 * d;
        return true;
    }
    if (t2 >= 0.f && t2 <= 1.f) {
        t_out = t2;
        intersection = o + t2 * d;
        return true;
    }

    return false;
}

std::vector<glm::vec2> generate_poisson_disk_grid(glm::vec2 center, float radius_max, float radius_min, int max_points = 10000)
{
    const float cell_size = radius_min / std::sqrt(2.f);
    const int grid_dim = int(std::ceil(2.f * radius_max / cell_size));
    std::vector<std::vector<int>> grid(grid_dim * grid_dim, std::vector<int>{});

    auto grid_index = [&](glm::vec2 p) -> glm::ivec2 {
        glm::vec2 offset = p - (center - glm::vec2(radius_max));
        return glm::clamp(glm::ivec2(offset / cell_size), glm::ivec2(0), glm::ivec2(grid_dim - 1));
    };

    std::vector<glm::vec2> points;
    std::vector<int> active;

    // 1. point initial au centre
    glm::vec2 first = center;
    points.push_back(first);
    active.push_back(0);
    glm::ivec2 gi = grid_index(first);
    grid[gi.x + gi.y * grid_dim].push_back(0);

    const int k = 100;

    while (!active.empty() && points.size() < max_points) {
        int idx = utils::rand(0, int(active.size()));
        glm::vec2 origin = points[active[idx]];
        bool found = false;

        for (int i = 0; i < k; ++i) {
            float angle = utils::rand(0.f, 2.f * 3.14159f);
            float dist = utils::rand(radius_min, 2.f * radius_min);
            glm::vec2 offset = dist * glm::vec2(std::cos(angle), std::sin(angle));
            glm::vec2 candidate = origin + offset;

            if (glm::distance(candidate, center) > radius_max)
                continue;

            glm::ivec2 cg = grid_index(candidate);
            bool valid = true;

            for (int y = -2; y <= 2 && valid; ++y) {
                for (int x = -2; x <= 2 && valid; ++x) {
                    glm::ivec2 neighbor = cg + glm::ivec2(x, y);
                    if (neighbor.x < 0 || neighbor.x >= grid_dim || neighbor.y < 0 || neighbor.y >= grid_dim)
                        continue;

                    for (int j : grid[neighbor.x + neighbor.y * grid_dim]) {
                        if (glm::distance(points[j], candidate) < radius_min) {
                            valid = false;
                            break;
                        }
                    }
                }
            }

            if (valid) {
                int new_index = int(points.size());
                points.push_back(candidate);
                active.push_back(new_index);
                grid[cg.x + cg.y * grid_dim].push_back(new_index);
                found = true;
                break;
            }
        }

        if (!found) {
            active[idx] = active.back();
            active.pop_back();
        }
    }

    std::cout << "Points generes avec poisson : " << points.size() << " / " << max_points << std::endl;
    return points;
}

struct Particle {
    glm::vec2 position;
    glm::vec2 velocity;
    float mass;
    float age = 0.f;
    float lifetime = 0.f;
    glm::vec4 color_start;
    glm::vec4 color_end;

    Particle()
    {
        glm::vec2 center = {0.0f, 0.0f};
        float radius_max = 0.3f;

        float boxAngle = utils::rand(0.f, 2.f * 3.14159f);
        float radius = radius_max * std::sqrt(utils::rand(0.f, 1.f));

        position = center + radius * glm::vec2(std::cos(boxAngle), std::sin(boxAngle));

        float angle = utils::rand(0.f, 2.f * 3.14159f);
        float speed = 0; //utils::rand(0.2f, 0.6f);
        mass = utils::rand(0.5f, 2.f);

        velocity = glm::vec2(std::cos(angle), std::sin(angle)) * speed;

        lifetime = utils::rand(5.f, 10.f);

        color_start = glm::vec4(utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), 1.f);
        color_end   = glm::vec4(utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), 1.f);
    }
};

void draw_parametric(std::function<glm::vec2(float)> const& parametric, int resolution = 100)
{
    for (int i = 0; i < resolution - 1; ++i) {
        float t1 = float(i) / (resolution - 1);
        float t2 = float(i + 1) / (resolution - 1);
        glm::vec2 p1 = parametric(t1);
        glm::vec2 p2 = parametric(t2);

        utils::draw_line(p1, p2, 0.008f, {1.f, 1.f, 1.f, 1.f});
    }
}

glm::vec2 bezier1_bernstein(glm::vec2 p0, glm::vec2 p1, float t)
{
    return (1 - t) * p0 + t * p1;
}
glm::vec2 bezier2_bernstein(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2, float t)
{
    float u = 1 - t;
    return u * u * p0 + 2 * u * t * p1 + t * t * p2;
}
glm::vec2 bezier3_bernstein(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, float t)
{
    float u = 1 - t;
    return u * u * u * p0
         + 3 * u * u * t * p1
         + 3 * u * t * t * p2
         + t * t * t * p3;
}

glm::vec2 bezier1_casteljau(glm::vec2 p0, glm::vec2 p1, float t)
{
    return (1 - t) * p0 + t * p1;
}
glm::vec2 bezier2_casteljau(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2, float t)
{
    glm::vec2 a = glm::mix(p0, p1, t);
    glm::vec2 b = glm::mix(p1, p2, t);
    return glm::mix(a, b, t);
}
glm::vec2 bezier3_casteljau(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, float t)
{
    glm::vec2 a = glm::mix(p0, p1, t);
    glm::vec2 b = glm::mix(p1, p2, t);
    glm::vec2 c = glm::mix(p2, p3, t);
    glm::vec2 d = glm::mix(a, b, t);
    glm::vec2 e = glm::mix(b, c, t);
    return glm::mix(d, e, t);
}

float find_nearest_t(std::function<glm::vec2(float)> bezier, glm::vec2 target)
{
    float best_t = 0.f;
    float best_dist2 = std::numeric_limits<float>::max();
    const int samples = 100;

    for (int i = 0; i <= samples; ++i) {
        float t_candidate = float(i) / samples;
        glm::vec2 diff = bezier(t_candidate) - target;
        float dist2 = glm::dot(diff, diff);
        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best_t = t_candidate;
        }
    }

    float t = best_t;
    const float step = 0.01f;
    const float epsilon = 1e-4f;
    const int max_iter = 100;

    for (int i = 0; i < max_iter; ++i) {
        float t1 = glm::clamp(t - epsilon, 0.f, 1.f);
        float t2 = glm::clamp(t + epsilon, 0.f, 1.f);

        glm::vec2 d1 = bezier(t1) - target;
        glm::vec2 d2 = bezier(t2) - target;

        float f1 = glm::dot(d1, d1);
        float f2 = glm::dot(d2, d2);

        float grad = (f2 - f1) / (2.f * epsilon);
        t = glm::clamp(t - step * grad, 0.f, 1.f);
    }

    return t;
}

float find_nearest_t_newton_safe(std::function<glm::vec2(float)> bezier, glm::vec2 target)
{
    const int steps = 50;
    float best_t = 0.f;
    float best_dist2 = std::numeric_limits<float>::max();

    for (int i = 0; i <= steps; ++i) {
        float t = float(i) / steps;
        glm::vec2 diff = bezier(t) - target;
        float dist2 = glm::dot(diff, diff);
        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best_t = t;
        }
    }

    float t = best_t;
    const int max_iter = 20;
    const float epsilon = 1e-5f;
    const float h = 1e-3f;

    for (int i = 0; i < max_iter; ++i) {
        glm::vec2 p = bezier(t);
        glm::vec2 diff = p - target;

        float t_forward = glm::clamp(t + h, 0.f, 1.f);
        float t_backward = glm::clamp(t - h, 0.f, 1.f);

        glm::vec2 p_plus = bezier(t_forward);
        glm::vec2 p_minus = bezier(t_backward);

        glm::vec2 d1 = (p_plus - p_minus) / (2.f * h);
        glm::vec2 d2 = (p_plus - 2.f * p + p_minus) / (h * h);

        float f_prime = 2.f * glm::dot(diff, d1);
        float f_double_prime = 2.f * (glm::dot(d1, d1) + glm::dot(diff, d2));

        if (std::abs(f_double_prime) < epsilon)
            break;

        float delta = f_prime / f_double_prime;
        t -= delta;
        t = glm::clamp(t, 0.f, 1.f);

        if (std::abs(delta) < epsilon)
            break;
    }

    return t;
}


int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    std::vector<Particle> particles;

    glm::vec2 p0 = {-0.8f,  0.0f};
    glm::vec2 p1 = {-0.4f,  0.6f};
    glm::vec2 p2 = { 0.4f, -0.6f};
    glm::vec2 p3 = { 0.8f,  0.0f};

    const int NUM_PARTICLES = 100;

    for (int i = 0; i < NUM_PARTICLES; ++i) {
        float t = float(i) / (NUM_PARTICLES - 1);
        glm::vec2 pos = bezier3_casteljau(p0, p1, p2, p3, t);

        float delta = 0.001f;
        glm::vec2 before = bezier3_casteljau(p0, p1, p2, p3, glm::clamp(t - delta, 0.f, 1.f));
        glm::vec2 after  = bezier3_casteljau(p0, p1, p2, p3, glm::clamp(t + delta, 0.f, 1.f));
        glm::vec2 tangent = glm::normalize(after - before);
        glm::vec2 normal = glm::vec2(-tangent.y, tangent.x);

        Particle p;
        p.position = pos;
        p.velocity = normal * 0.3f;
        p.mass = 1.0f;
        p.lifetime = 10.0f;
        p.color_start = {1.f, 0.f, 0.f, 1.f};
        p.color_end   = {1.f, 0.f, 0.f, 1.f};

        particles.push_back(p);
    }

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        float dt = gl::delta_time_in_seconds();

        draw_parametric([&](float t) {
            return bezier3_casteljau(p0, p1, p2, p3, t);
        });

        glm::vec2 mouse = (gl::mouse_position());
        utils::draw_disk(mouse, 0.015f, {1.f, 1.f, 0.f, 1.f}); 
        

        auto bezier = [&](float t) {
            return bezier3_casteljau(p0, p1, p2, p3, t);
        };

        //float nearest_t = find_nearest_t(bezier, mouse);
        float nearest_t = find_nearest_t_newton_safe(bezier, mouse);

        glm::vec2 closest = bezier(nearest_t);
        utils::draw_disk(closest, 0.02f, {0.f, 1.f, 0.f, 1.f}); // cercle vert au point le + proche
        utils::draw_line(mouse, closest, 0.002f, {0.f, 1.f, 1.f, 1.f});

        for (auto& p : particles) {
            glm::vec2 total_force{};
            total_force += glm::vec2(0.f, -1.0f * p.mass);  // gravité vers le bas

            glm::vec2 gravity = glm::vec2(0.f, -1.00f * p.mass);

            auto bezier = [&](float t) {
                return bezier3_casteljau(p0, p1, p2, p3, t);
            };

            float t_near = find_nearest_t(bezier, p.position);
            glm::vec2 nearest_point = bezier(t_near);

            float delta = 0.001f;
            glm::vec2 before = bezier(glm::clamp(t_near - delta, 0.f, 1.f));
            glm::vec2 after  = bezier(glm::clamp(t_near + delta, 0.f, 1.f));
            glm::vec2 tangent = glm::normalize(after - before);
            glm::vec2 normal = glm::vec2(-tangent.y, tangent.x);

            glm::vec2 dir = p.position - nearest_point;
            float dist = glm::length(dir);

            if (dist < 0.2f && dist > 0.0001f) {
                glm::vec2 dir_norm = glm::normalize(dir);
                float strength = (1.f - dist / 0.2f);  // décroissant
                glm::vec2 force = dir_norm * strength * 5.f;  // **repousse**
                total_force += force;
            }


            /*
            float k = 1.0f;
            //glm::vec2 air_resistance = -k * p.velocity;
           // total_force += air_resistance;
           */

            glm::vec2 acceleration = total_force / p.mass;
            p.velocity += acceleration * dt;

            p.position += p.velocity * dt;
            p.age += dt;
        }

        for (const auto& p : particles) {
            /*
            float shrink_duration = 2.f;
            float time_left = p.lifetime - p.age;
            float shrink_factor = glm::clamp(time_left / shrink_duration, 0.f, 1.f);

            float pulse = glm::sin(p.age * 20.f) * 0.5f + 0.5f;
            pulse = glm::mix(1.f, pulse, 1.f - shrink_factor);

            float radius = 0.02f * shrink_factor * pulse;
            */
            float radius = 0.01f; // taille fixe sans shrink

            float t = glm::clamp(p.age / p.lifetime, 0.f, 1.f);
            glm::vec4 color = glm::mix(p.color_start, p.color_end, t);
            // color.a *= shrink_factor; // transparence désactivée

            utils::draw_disk(p.position, radius, color);
        }

        /*
        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                [](const Particle& p) {
                    return p.age > p.lifetime;
                }),
            particles.end()
        );
        */
    }
}
