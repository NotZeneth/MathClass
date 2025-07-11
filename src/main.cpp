#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

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

    // On résout||o + d·t - c||² = r²
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
        
        glm::vec2 spawn_min = { 0.2f, -0.8f};
        glm::vec2 spawn_max = { 0.8f,  0.8f};

        position = {
            utils::rand(spawn_min.x, spawn_max.x),
            utils::rand(spawn_min.y, spawn_max.y)
        };

        float angle = utils::rand(0.f, 2.f * 3.14159f);
        float speed = 0; //utils::rand(0.2f, 0.6f);
        mass = utils::rand(0.5f, 2.f);

        velocity = glm::vec2(std::cos(angle), std::sin(angle)) * speed;

        lifetime = utils::rand(5.f, 10.f);

        color_start = glm::vec4(utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), 1.f);
        color_end   = glm::vec4(utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), utils::rand(0.2f, 1.f), 1.f);
    }
};

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    std::vector<Particle> particles;
    for (int i = 0; i < 100; ++i) {
        particles.emplace_back();
    }

    std::vector<Segment> walls = {
        Segment(glm::vec2(-1.f, -0.8f), glm::vec2(1.f, -0.5f)),
        Segment(glm::vec2(-0.5f, -1.f), glm::vec2(-0.2f, 1.f)),
        Segment(glm::vec2(-1.f, 0.8f), glm::vec2(1.f, 0.9f)),
        Segment(glm::vec2(-gl::window_aspect_ratio(), -1.f), glm::vec2(gl::window_aspect_ratio(), -1.f)),
        Segment(glm::vec2(-gl::window_aspect_ratio(),  1.f), glm::vec2(gl::window_aspect_ratio(),  1.f)),
        Segment(glm::vec2(-gl::window_aspect_ratio(), -1.f), glm::vec2(-gl::window_aspect_ratio(), 1.f)),
        Segment(glm::vec2(gl::window_aspect_ratio(), -1.f), glm::vec2(gl::window_aspect_ratio(), 1.f)),
    };

    std::vector<Circle> circles = {
        Circle(glm::vec2(0.f, 0.f), 0.3f),
        Circle(glm::vec2(0.5f, -0.3f), 0.2f)
    };

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        float dt = gl::delta_time_in_seconds();

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

            for (const auto& wall : walls) {
                glm::vec2 intersection;
                float t;
                if (intersect_segments(old_position, new_position, wall.a, wall.b, intersection, t)) {
                    glm::vec2 wall_dir = glm::normalize(wall.b - wall.a);
                    glm::vec2 normal = glm::vec2(-wall_dir.y, wall_dir.x);

                    p.velocity = glm::reflect(p.velocity, normal);

                    float distance_left = glm::length(new_position - intersection);
                    p.position = intersection + p.velocity * (distance_left / glm::length(p.velocity));

                    collided = true;
                    break;
                }
            }

            for (const auto& circle : circles) {
                glm::vec2 intersection;
                float t;

                if (intersect_segment_circle(old_position, new_position, circle.center, circle.radius, intersection, t)) {
                    glm::vec2 normal = glm::normalize(intersection - circle.center);
                    p.velocity = glm::reflect(p.velocity, normal);

                    float distance_left = glm::length(new_position - intersection);
                    p.position = intersection + p.velocity * (distance_left / glm::length(p.velocity));

                    collided = true;
                    break;
                }
            }

            if (!collided)
                p.position = new_position;

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
            float radius = 0.02f; // taille fixe sans shrink

            float t = glm::clamp(p.age / p.lifetime, 0.f, 1.f);
            glm::vec4 color = glm::mix(p.color_start, p.color_end, t);
            // color.a *= shrink_factor; // transparence désactivée

            utils::draw_disk(p.position, radius, color);
        }

        for (const auto& wall : walls) {
            utils::draw_line(wall.a, wall.b, 0.01f, {1.f, 1.f, 1.f, 1.f});
        }

        for (const auto& circle : circles) {
            utils::draw_disk(circle.center, circle.radius, {1.f, 0.f, 0.f, 1.f});
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
