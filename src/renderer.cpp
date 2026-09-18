#include "renderer.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

// ============================================================================
// Shaders
// ============================================================================

// One shader body for both targets, with the dialect supplied as a separate
// source string.  The bodies below are already legal in both GLSL 3.30 and
// GLSL ES 3.00 — `layout(location = ...)` on vertex inputs and a user-declared
// fragment output are in both — so only the version line and the ES precision
// qualifiers differ, and they are the whole of the difference.
//
// ES has no default precision for `float` in a fragment shader: omit the
// qualifier and the shader fails to compile.  Under WebGL that failure is
// silent, and presents as a black canvas with nothing but a console message.
// The vertex stage does default to highp and would compile without it, but the
// prologue is one string for both stages — stating the precision that stage
// already has is free, and splitting the prologue in two to avoid saying it
// would only invite the halves to drift.
#ifdef __EMSCRIPTEN__
static const char* shader_prologue =
    "#version 300 es\n"
    "precision highp float;\n";
#else
static const char* shader_prologue = "#version 330 core\n";
#endif

static const char* vertex_shader_src = R"(
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec4 aColor;
uniform mat4 uVP;
out vec4 vColor;
void main() {
    gl_Position = uVP * vec4(aPos, 1.0);
    vColor = aColor;
}
)";

static const char* fragment_shader_src = R"(
in vec4 vColor;
out vec4 FragColor;
void main() {
    FragColor = vColor;
}
)";

// ============================================================================
// LineRenderer
// ============================================================================

GLuint LineRenderer::compile_shader(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    const char* sources[] = {shader_prologue, source};
    glShaderSource(shader, 2, sources, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        std::fprintf(stderr, "Shader compilation failed: %s\n", log);
        std::exit(1);
    }
    return shader;
}

GLuint LineRenderer::link_program(GLuint vs, GLuint fs) {
    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char log[512];
        glGetProgramInfoLog(program, 512, nullptr, log);
        std::fprintf(stderr, "Program link failed: %s\n", log);
        std::exit(1);
    }

    glDeleteShader(vs);
    glDeleteShader(fs);
    return program;
}

LineRenderer::LineRenderer() {
    GLuint vs = compile_shader(GL_VERTEX_SHADER, vertex_shader_src);
    GLuint fs = compile_shader(GL_FRAGMENT_SHADER, fragment_shader_src);
    shader_program_ = link_program(vs, fs);
    vp_loc_ = glGetUniformLocation(shader_program_, "uVP");

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);

    // Position (3 floats)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ColorVertex),
                          (void*)offsetof(ColorVertex, x));
    glEnableVertexAttribArray(0);

    // Color (4 floats)
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(ColorVertex),
                          (void*)offsetof(ColorVertex, r));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}

LineRenderer::~LineRenderer() {
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (shader_program_) glDeleteProgram(shader_program_);
}

void LineRenderer::set_vp(const Eigen::Matrix4f& vp) {
    glUseProgram(shader_program_);
    glUniformMatrix4fv(vp_loc_, 1, GL_FALSE, vp.data());
}

void LineRenderer::begin() {
    batches_.clear();
    current_batch_ = nullptr;
    current_width_ = -1.0f;
    tri_vertices_.clear();
}

void LineRenderer::ensure_batch(float width) {
    if (current_batch_ && std::abs(current_width_ - width) < 0.01f) return;
    batches_.push_back({{}, width});
    current_batch_ = &batches_.back();
    current_width_ = width;
}

void LineRenderer::line(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                        const std::array<float, 4>& c, float width) {
    ensure_batch(width);
    current_batch_->vertices.push_back(
        {(float)a.x(), (float)a.y(), (float)a.z(), c[0], c[1], c[2], c[3]});
    current_batch_->vertices.push_back(
        {(float)b.x(), (float)b.y(), (float)b.z(), c[0], c[1], c[2], c[3]});
}

void LineRenderer::circle(const Eigen::Vector3d& center, double radius,
                          const Eigen::Vector3d& normal,
                          const std::array<float, 4>& color,
                          int segments, float width) {
    // Build orthonormal basis on the plane
    Eigen::Vector3d n = normal.normalized();
    Eigen::Vector3d u, v;
    if (std::abs(n.x()) < 0.9) {
        u = Eigen::Vector3d::UnitX().cross(n).normalized();
    } else {
        u = Eigen::Vector3d::UnitY().cross(n).normalized();
    }
    v = n.cross(u);

    for (int i = 0; i < segments; ++i) {
        double t0 = 2.0 * M_PI * i / segments;
        double t1 = 2.0 * M_PI * (i + 1) / segments;
        Eigen::Vector3d p0 = center + radius * (std::cos(t0) * u + std::sin(t0) * v);
        Eigen::Vector3d p1 = center + radius * (std::cos(t1) * u + std::sin(t1) * v);
        line(p0, p1, color, width);
    }
}

void LineRenderer::point(const Eigen::Vector3d& pos,
                         const std::array<float, 4>& color,
                         float size, float width) {
    Eigen::Vector3d dx(size, 0, 0), dy(0, size, 0), dz(0, 0, size);
    line(pos - dx, pos + dx, color, width);
    line(pos - dy, pos + dy, color, width);
    line(pos - dz, pos + dz, color, width);
}

void LineRenderer::axes(const Eigen::Vector3d& origin, double length) {
    line(origin, origin + Eigen::Vector3d(length, 0, 0), {1, 0, 0, 1}, 2.0f);
    line(origin, origin + Eigen::Vector3d(0, length, 0), {0, 1, 0, 1}, 2.0f);
    line(origin, origin + Eigen::Vector3d(0, 0, length), {0, 0, 1, 1}, 2.0f);
}

void LineRenderer::disc(const Eigen::Vector3d& center, double radius,
                        const Eigen::Vector3d& normal,
                        const std::array<float, 4>& color, int segments) {
    Eigen::Vector3d n = normal.normalized();
    Eigen::Vector3d u, v;
    if (std::abs(n.x()) < 0.9) {
        u = Eigen::Vector3d::UnitX().cross(n).normalized();
    } else {
        u = Eigen::Vector3d::UnitY().cross(n).normalized();
    }
    v = n.cross(u);

    ColorVertex cv_center = {(float)center.x(), (float)center.y(), (float)center.z(),
                             color[0], color[1], color[2], color[3]};

    for (int i = 0; i < segments; ++i) {
        double t0 = 2.0 * M_PI * i / segments;
        double t1 = 2.0 * M_PI * (i + 1) / segments;
        Eigen::Vector3d p0 = center + radius * (std::cos(t0) * u + std::sin(t0) * v);
        Eigen::Vector3d p1 = center + radius * (std::cos(t1) * u + std::sin(t1) * v);

        tri_vertices_.push_back(cv_center);
        tri_vertices_.push_back({(float)p0.x(), (float)p0.y(), (float)p0.z(),
                                 color[0], color[1], color[2], color[3]});
        tri_vertices_.push_back({(float)p1.x(), (float)p1.y(), (float)p1.z(),
                                 color[0], color[1], color[2], color[3]});
    }
}

void LineRenderer::flush() {
    glUseProgram(shader_program_);
    glBindVertexArray(vao_);

    // Draw filled triangles with depth-write disabled (see-through surfaces)
    if (!tri_vertices_.empty()) {
        glDepthMask(GL_FALSE);  // Don't write to depth buffer
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glBufferData(GL_ARRAY_BUFFER,
                     tri_vertices_.size() * sizeof(ColorVertex),
                     tri_vertices_.data(), GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(tri_vertices_.size()));
        glDepthMask(GL_TRUE);   // Re-enable for lines
    }

    // Draw lines on top.
    //
    // Line width under WebGL2 is the driver's to honour: issue #11 assumed it
    // was silently ignored above 1.0, and on this machine's Chrome it is not —
    // ALIASED_LINE_WIDTH_RANGE reads [1, 10] and the widths below rasterize.
    // Elsewhere the range is [1, 1] and every line comes out one pixel.  Both
    // are conforming, so width is decorative here and carries no information:
    // every element of the scene is told apart by colour, and reads correctly
    // at a uniform 1 px.  Measured on the web build for #15.
    for (const auto& batch : batches_) {
        if (batch.vertices.empty()) continue;

        glLineWidth(batch.width);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glBufferData(GL_ARRAY_BUFFER,
                     batch.vertices.size() * sizeof(ColorVertex),
                     batch.vertices.data(), GL_DYNAMIC_DRAW);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(batch.vertices.size()));
    }

    glBindVertexArray(0);
}

// ============================================================================
// Camera
// ============================================================================

Eigen::Vector3f OrbitCamera::eye_position() const {
    float az = azimuth * M_PI / 180.0f;
    float el = elevation * M_PI / 180.0f;
    float cx = distance * std::cos(el) * std::cos(az);
    float cy = distance * std::cos(el) * std::sin(az);
    float cz = distance * std::sin(el);
    return target + Eigen::Vector3f(cx, cy, cz);
}

Eigen::Matrix4f OrbitCamera::view_matrix() const {
    return look_at(eye_position(), target, Eigen::Vector3f(0, 0, 1));
}

// ============================================================================
// Matrix Helpers
// ============================================================================

Eigen::Matrix4f perspective(float fov_deg, float aspect, float near, float far) {
    float fov = fov_deg * M_PI / 180.0f;
    float t = std::tan(fov / 2.0f);

    Eigen::Matrix4f m = Eigen::Matrix4f::Zero();
    m(0, 0) = 1.0f / (aspect * t);
    m(1, 1) = 1.0f / t;
    m(2, 2) = -(far + near) / (far - near);
    m(2, 3) = -2.0f * far * near / (far - near);
    m(3, 2) = -1.0f;
    return m;
}

Eigen::Matrix4f look_at(const Eigen::Vector3f& eye, const Eigen::Vector3f& target,
                         const Eigen::Vector3f& up) {
    Eigen::Vector3f f = (target - eye).normalized();
    Eigen::Vector3f s = f.cross(up).normalized();
    Eigen::Vector3f u = s.cross(f);

    Eigen::Matrix4f m = Eigen::Matrix4f::Identity();
    m(0, 0) = s.x(); m(0, 1) = s.y(); m(0, 2) = s.z(); m(0, 3) = -s.dot(eye);
    m(1, 0) = u.x(); m(1, 1) = u.y(); m(1, 2) = u.z(); m(1, 3) = -u.dot(eye);
    m(2, 0) = -f.x(); m(2, 1) = -f.y(); m(2, 2) = -f.z(); m(2, 3) = f.dot(eye);
    return m;
}
