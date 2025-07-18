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

        utils::draw_line(p1, p2, 0.002f, {1.f, 1.f, 1.f, 1.f});
    }
}

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    std::vector<Particle> particles;

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        float dt = gl::delta_time_in_seconds();

        // Cercle
        draw_parametric([](float t) {
            float angle = t * 2.f * 3.14159f;
            return glm::vec2(std::cos(angle), std::sin(angle)) * 0.4f;
        });

        // Coeur
        draw_parametric([](float t) {
            float angle = t * 2.f * 3.14159f;
            float x = 0.16f * std::sin(angle) * std::sin(angle) * std::sin(angle);
            float y = 0.13f * std::cos(angle) - 0.05f * std::cos(2 * angle)
                    - 0.02f * std::cos(3 * angle) - 0.01f * std::cos(4 * angle);
            return glm::vec2(x, y);
        });

        for (auto& p : particles) {
            glm::vec2 total_force{};

            glm::vec2 gravity = glm::vec2(0.f, -1.00f * p.mass);

            /*
            float k = 1.0f;
            //glm::vec2 air_resistance = -k * p.velocity;
           // total_force += air_resistance;
           */

            glm::vec2 acceleration = total_force / p.mass;
            p.velocity += acceleration * dt;

            glm::vec2 old_position = p.position;
            glm::vec2 new_position = p.position + p.velocity * dt;

            bool collided = false;

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
            float radius = 0.001f; // taille fixe sans shrink

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
