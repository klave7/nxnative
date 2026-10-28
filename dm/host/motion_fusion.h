// Controller motion (gyro + accelerometer) for the six-axis sensor (hle_hid.c). Pure C, no host state, so tests/t_motion.c checks it on any machine.
// Input: what the front end reads from GameController (GCMotion), in the controller's own frame as Apple reports it: acceleration in G (gravity direction, as CoreMotion:
// (0, 0, -1) lying face up), rotation rate in rad/s. Output: the Switch's SixAxisSensorState fields: acceleration in G, angular velocity in revolutions/s, angle in revolutions,
// direction = the controller's three axes in the world (rows). Switch frame: x right, y forward (away from the player), z up out of the face; flat at rest acc = (0, 0, -1).
// The frame map (DM_MOTION_AXES, default "x,y,z") is a signed permutation from Apple's axes to the Switch's, applied to acceleration and rotation alike (rotation times its
// determinant, so a mirror stays a proper rotation); DM_MOTION_ACC_SIGN=-1 flips only the acceleration, for a source that reports the specific force (+1 G up at rest).
// Orientation: the gyro is integrated (R = R exp(w dt)), and tilt drift is pulled back toward the measured gravity (Mahony, gain k) while |acc| is near 1 G; yaw is free to drift,
// as on the console without a recenter.
#pragma once
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct { int src[3]; float sign[3]; float acc_sign; float det; } motion_map;
typedef struct {
    float R[3][3];        // controller -> world (columns = the controller's axes in the world)
    float angle[3];       // integrated revolutions per axis
    float acc[3], gyro[3];// latest mapped sample: G, revolutions / s
    int primed;
} motion_state;

// "x,y,z" with optional '-' per entry, e.g. "x,-z,y". Returns 0 and the identity on a malformed string
static inline int motion_map_parse(motion_map *m, const char *s) {
    motion_map id = {{0, 1, 2}, {1, 1, 1}, 1, 1};
    *m = id;
    if (!s || !*s) return 1;
    motion_map r = id; int used = 0;
    for (int i = 0; i < 3; i++) {
        while (*s == ' ') s++;
        float sg = 1; if (*s == '-') { sg = -1; s++; } else if (*s == '+') s++;
        int a = *s == 'x' || *s == 'X' ? 0 : *s == 'y' || *s == 'Y' ? 1 : *s == 'z' || *s == 'Z' ? 2 : -1;
        if (a < 0 || used >> a & 1) return 0;
        used |= 1 << a; r.src[i] = a; r.sign[i] = sg; s++;
        while (*s == ' ') s++;
        if (i < 2) { if (*s != ',') return 0; s++; }
    }
    if (*s) return 0;
    // determinant of the signed permutation: sign product times the permutation's parity
    int inv = (r.src[0] > r.src[1]) + (r.src[0] > r.src[2]) + (r.src[1] > r.src[2]);
    r.det = r.sign[0] * r.sign[1] * r.sign[2] * (inv & 1 ? -1.f : 1.f);
    *m = r;
    return 1;
}
static inline void motion_map_apply(const motion_map *m, const float in[3], float out[3], float scale) {
    for (int i = 0; i < 3; i++) out[i] = m->sign[i] * in[m->src[i]] * scale;
}
static inline void motion_reset(motion_state *s) {
    memset(s, 0, sizeof *s);
    for (int i = 0; i < 3; i++) s->R[i][i] = 1;
    s->acc[2] = -1;
}
static inline void m3_mul(float o[3][3], const float a[3][3], const float b[3][3]) {
    float t[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) t[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    memcpy(o, t, sizeof t);
}
static inline void m3_orthonormalize(float R[3][3]) {   // Gram-Schmidt on the columns
    float c[3][3]; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) c[j][i] = R[i][j];
    for (int k = 0; k < 3; k++) {
        for (int p = 0; p < k; p++) { float d = c[k][0] * c[p][0] + c[k][1] * c[p][1] + c[k][2] * c[p][2]; for (int i = 0; i < 3; i++) c[k][i] -= d * c[p][i]; }
        float n = sqrtf(c[k][0] * c[k][0] + c[k][1] * c[k][1] + c[k][2] * c[k][2]); if (n < 1e-6f) n = 1;
        for (int i = 0; i < 3; i++) c[k][i] /= n;
    }
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) R[i][j] = c[j][i];
}
// One sample: raw Apple acceleration (G) and rotation rate (rad/s), dt in seconds since the previous one (0 for the first). k = tilt correction gain (rad/s per unit error)
static inline void motion_update(motion_state *s, const motion_map *m, const float acc_raw[3], const float rot_raw[3], float dt, float k) {
    float a[3], w[3];
    motion_map_apply(m, acc_raw, a, m->acc_sign);
    motion_map_apply(m, rot_raw, w, m->det);           // rad/s in the Switch frame
    if (!s->primed || dt <= 0 || dt > 0.25f) {         // first sample or a gap: no integration step
        if (!s->primed) {                              // start level with the measured gravity, facing forward
            float n = sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
            if (n > 0.5f && n < 1.5f) {
                // world down (0,0,-1) = R g_body  ->  R's third row = -g_body / |g|
                float z[3] = {-a[0] / n, -a[1] / n, -a[2] / n}, y0[3] = {0, 1, 0};
                // Upright phones have gravity along y: choose another heading reference before projecting it onto the level plane.
                if (fabsf(z[1]) > 0.99f) { y0[0] = 1; y0[1] = 0; }
                float d = y0[0] * z[0] + y0[1] * z[1] + y0[2] * z[2];
                float y[3] = {y0[0] - d * z[0], y0[1] - d * z[1], y0[2] - d * z[2]};
                float yn = sqrtf(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
                if (yn > 1e-3f) {
                    for (int i = 0; i < 3; i++) y[i] /= yn;
                    float x[3] = {y[1] * z[2] - y[2] * z[1], y[2] * z[0] - y[0] * z[2], y[0] * z[1] - y[1] * z[0]};
                    for (int j = 0; j < 3; j++) { s->R[0][j] = x[j]; s->R[1][j] = y[j]; s->R[2][j] = z[j]; }
                }
            }
            s->primed = 1;
        }
    } else {
        // tilt correction: predicted gravity in the body = R^T (0,0,-1) = -(third row of R)
        float n = sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if (k > 0 && n > 0.8f && n < 1.2f) {
            float gm[3] = {a[0] / n, a[1] / n, a[2] / n}, gp[3] = {-s->R[2][0], -s->R[2][1], -s->R[2][2]};
            float e[3] = {gm[1] * gp[2] - gm[2] * gp[1], gm[2] * gp[0] - gm[0] * gp[2], gm[0] * gp[1] - gm[1] * gp[0]};
            for (int i = 0; i < 3; i++) w[i] += k * e[i];
        }
        float th = sqrtf(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) * dt;
        if (th > 1e-9f) {
            float u[3] = {w[0] * dt / th, w[1] * dt / th, w[2] * dt / th}, c = cosf(th), sn = sinf(th), v = 1 - c;
            float E[3][3] = {
                {c + u[0] * u[0] * v,        u[0] * u[1] * v - u[2] * sn, u[0] * u[2] * v + u[1] * sn},
                {u[1] * u[0] * v + u[2] * sn, c + u[1] * u[1] * v,        u[1] * u[2] * v - u[0] * sn},
                {u[2] * u[0] * v - u[1] * sn, u[2] * u[1] * v + u[0] * sn, c + u[2] * u[2] * v}};
            m3_mul(s->R, s->R, E);
            m3_orthonormalize(s->R);
        }
    }
    // reported angular velocity: the gyro itself (revolutions / s), without the tilt correction term
    float wraw[3]; motion_map_apply(m, rot_raw, wraw, m->det / (2 * (float)M_PI));
    for (int i = 0; i < 3; i++) {
        s->gyro[i] = wraw[i]; s->acc[i] = a[i];
        if (dt > 0 && dt <= 0.25f) s->angle[i] += wraw[i] * dt;
    }
}
// The SixAxisSensorState fields from the state: acc[3], gyro[3], angle[3], dir[9] (rows = the controller's x, y, z axes in the world)
static inline void motion_fill(const motion_state *s, float acc[3], float gyro[3], float angle[3], float dir[9]) {
    for (int i = 0; i < 3; i++) { acc[i] = s->acc[i]; gyro[i] = s->gyro[i]; angle[i] = s->angle[i]; }
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) dir[3 * i + j] = s->R[j][i];
}
// The phone's own sensors (CoreMotion: x to the right of the screen in portrait, y to its top, z out of the screen) in the Switch frame for how the screen is turned, so a phone
// held like a Switch in handheld (or mounted in a controller) is the console: x to the right of the picture, y to its top, z toward the player. orient = UIInterfaceOrientation:
// 1 portrait, 2 upside down, 3 landscape right (home indicator right: the phone's top to the left), 4 landscape left (the phone's top to the right). Proper rotations only,
// so the same map serves acceleration and rotation rate
static inline void motion_phone_to_switch(int orient, const float in[3], float out[3]) {
    switch (orient) {
    case 3: out[0] = -in[1]; out[1] = in[0]; out[2] = in[2]; break;
    case 4: out[0] = in[1]; out[1] = -in[0]; out[2] = in[2]; break;
    case 2: out[0] = -in[0]; out[1] = -in[1]; out[2] = in[2]; break;
    default: out[0] = in[0]; out[1] = in[1]; out[2] = in[2]; break;
    }
}
// Which way the phone is really held, from gravity (CoreMotion's gravity vector: toward the ground, in g, in the phone's axes). The screen's orientation says how the picture is turned
// against the phone's body, which a phone whose screen changes under the player (a foldable opened or closed) can leave as it was while the body is held 90 degrees differently,
// so its gyro would turn the wrong way. Held upright, gravity runs down the picture (Switch -y) in the right orientation, whichever the screen says. The orientation that puts gravity
// most nearly down the picture is a candidate; once the same one has been the candidate (clearly: more than 0.6 g down the picture, so not lying flat) for `samples` in a row, it
// replaces the screen's orientation; when it is the screen's own again the replacement ends. On a phone whose screen follows its body the candidate is always the screen's, so nothing changes
typedef struct { int over, cand, run; } motion_orient;   // over: the orientation that replaces the screen's (0 none), cand: the current candidate, run: samples it has held
static inline int motion_orient_from_gravity(const float g[3], float *down) {
    int best = 1; float bd = 2;
    for (int o = 1; o <= 4; o++) { float s[3]; motion_phone_to_switch(o, g, s); if (s[1] < bd) { bd = s[1]; best = o; } }
    if (down) *down = bd;
    return best;
}
// One sample: returns the orientation to map the sensors with (1..4); iface is the screen's
static inline int motion_orient_update(motion_orient *t, int iface, const float g[3], int samples) {
    float down; int c = motion_orient_from_gravity(g, &down);
    if (down > -0.6f) t->run = 0;   // lying flat or in between: keep what was decided
    else {
        if (c != t->cand) { t->cand = c; t->run = 0; }
        if (++t->run >= samples) t->over = c == iface ? 0 : c;
    }
    return t->over ? t->over : iface;
}
