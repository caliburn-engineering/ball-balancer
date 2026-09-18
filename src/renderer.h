#pragma once

// GL headers: GLAD on the desktop, Emscripten's own ES 3.0 headers on the
// web.  Everything below is the portable subset WebGL2 also provides.
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <glad/glad.h>
#endif
#include <Eigen/Dense>
#include <vector>
#include <array>

/// Orbit camera around a target point
struct OrbitCamera {
    float azimuth   = 45.0f;    // Horizontal angle [deg]
    float elevation = 30.0f;    // Vertical angle [deg] (0 = level, 90 = top-down)
    float distance  = 0.8f;     // Distance from target [m]
    Eigen::Vector3f target{0.0f, 0.0f, 0.15f};  // Look-at point

    Eigen::Matrix4f view_matrix() const;
    Eigen::Vector3f eye_position() const;
};

/// Vertex with position and color
struct ColorVertex {
    float x, y, z;
    float r, g, b, a;
};

/// Simple line renderer over the GL subset shared by OpenGL 3.3 core and
/// WebGL2: vertex arrays, buffer objects and array draws.
class LineRenderer {
public:
    LineRenderer();
    ~LineRenderer();

    // Non-copyable
    LineRenderer(const LineRenderer&) = delete;
    LineRenderer& operator=(const LineRenderer&) = delete;

    /// Set the view-projection matrix (call once per frame)
    void set_vp(const Eigen::Matrix4f& vp);

    /// Begin a new frame of line drawing
    void begin();

    /// Add a line segment
    void line(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
              const std::array<float, 4>& color, float width = 2.0f);

    /// Add a circle (approximated as line segments)
    void circle(const Eigen::Vector3d& center, double radius,
                const Eigen::Vector3d& normal, const std::array<float, 4>& color,
                int segments = 48, float width = 1.5f);

    /// Add a point (drawn as a small cross)
    void point(const Eigen::Vector3d& pos, const std::array<float, 4>& color,
               float size = 0.005f, float width = 3.0f);

    /// Draw a coordinate frame (RGB = XYZ)
    void axes(const Eigen::Vector3d& origin, double length = 0.05);

    /// Add a filled disc (triangle fan with transparency).
    /// Rendered with depth-write disabled so geometry behind remains visible.
    void disc(const Eigen::Vector3d& center, double radius,
              const Eigen::Vector3d& normal, const std::array<float, 4>& color,
              int segments = 48);

    /// Flush all accumulated geometry to GPU
    void flush();

private:
    GLuint shader_program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint  vp_loc_ = -1;

    struct LineBatch {
        std::vector<ColorVertex> vertices;
        float width;
    };
    std::vector<LineBatch> batches_;
    LineBatch* current_batch_ = nullptr;
    float current_width_ = 2.0f;

    // Triangle batches for filled geometry
    std::vector<ColorVertex> tri_vertices_;

    void ensure_batch(float width);
    static GLuint compile_shader(GLenum type, const char* source);
    static GLuint link_program(GLuint vs, GLuint fs);
};

/// Build a perspective projection matrix
Eigen::Matrix4f perspective(float fov_deg, float aspect, float near, float far);

/// Build a look-at view matrix
Eigen::Matrix4f look_at(const Eigen::Vector3f& eye, const Eigen::Vector3f& target,
                         const Eigen::Vector3f& up);
