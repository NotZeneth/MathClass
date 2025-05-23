#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

struct Particle {
    glm::vec2 position;
    glm::vec2 velocity;
    float mass;

    Particle()
    {
        position = {
            utils::rand(-gl::window_aspect_ratio(), gl::window_aspect_ratio()),
            utils::rand(-1.f, 1.f)
        };

        // random speed nd stats
        float angle = utils::rand(0.f, 2.f * 3.14159f); //environ pi
        float speed = utils::rand(0.2f, 0.6f);
        mass = utils::rand(0.5f, 2.f); // rand between 0.5 et 2

        velocity = glm::vec2(std::cos(angle), std::sin(angle)) * speed;
    }
};


int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // Vect de particule
    std::vector<Particle> particles;

    for (int i = 0; i < 100; ++i) {
        particles.emplace_back();
    }

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        float dt = gl::delta_time_in_seconds();

        // Update de la pos
        for (auto& p : particles) {
            
            glm::vec2 total_force{};

            //Gravity
            glm::vec2 gravity = glm::vec2(0.f, -1.00f * p.mass); 
            //total_force += gravity;

            float k = 1.0f; // Resistance
            glm::vec2 air_resistance = -k * p.velocity;
            total_force += air_resistance;

            //Acceleration
            glm::vec2 acceleration = total_force / p.mass;
            p.velocity += acceleration * dt;

            p.position += p.velocity * dt; // Update finale de la pos
            
        }

        // Render des particules
        for (const auto& p : particles) {
            utils::draw_disk(p.position, 0.02f, glm::vec4(1.f, 1.f, 1.f, 1.f));
        }

    }
}