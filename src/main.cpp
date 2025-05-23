#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

struct Particle {
    glm::vec2 position;

    Particle()
        : position{utils::rand(-gl::window_aspect_ratio(), gl::window_aspect_ratio()),
                   utils::rand(-1.f, 1.f)} {}
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

        // Render des particules
        for (const auto& p : particles) {
            utils::draw_disk(p.position, 0.02f, glm::vec4(1.f, 1.f, 1.f, 1.f));
        }

    }
}