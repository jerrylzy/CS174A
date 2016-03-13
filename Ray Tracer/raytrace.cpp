//
//  raytrace.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

//  Size of the Picture

int g_width;
int g_height;

//  Global Constants for Intersection Function
const float RECOLLISION = 0.0001f;
const float MIN_DIST = 1.0f;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

//  Structure of Properties for Sphere

struct Sphere
{
    vec4   pos;
    vec4   color;
    mat4   trans;
    mat4   inverseTrans;
    float  ambient;
    float  diffuse;
    float  specular;
    float  reflection;
    float  specExp;
    string name;
};

//  Structure of Properties for Light

struct Light
{
    vec4   pos;
    vec4   intensity;
    string name;
};

//  Loaded Properties

vector<vec4>    g_colors;   //  global color array
vector<Sphere>  g_sphere;   //  container for spheres
vector<Light>   g_light;    //  container for lights

//  Planes
float  g_left;
float  g_right;
float  g_top;
float  g_bottom;
float  g_near;

//  Global color constants
vec4   g_back;
vec4   g_ambient;

//  File name
string filename;


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

inline
float normalizeCondition(float value)
{
    return value > 1.0f ? 1.0f : value;
}

//  Normalize colors
inline
vec4 normalizeColorComponent(const vec4& value)
{
    return vec4(normalizeCondition(value.x), normalizeCondition(value.y),
                normalizeCondition(value.z), 1.0f);
}

//  Parse Line Utilities

const int NUM_LABELS  = 11;
const string LABELS[] = { "NEAR", "LEFT", "RIGHT", "BOTTOM", "TOP", "SPHERE",
    "LIGHT", "BACK", "AMBIENT", "OUTPUT", "RES" };

//  Return the first 3 components of a vec4
inline
vec3 toVec3(vec4 in)
{
    return vec3(in.x, in.y, in.z);
}

//  Get the label for a specific string
inline
size_t getLabel(const string& s)
{
    return find(LABELS, LABELS + NUM_LABELS, s) - LABELS;
}

// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

//  Load data from each line of strings

void parseLine(const vector<string>& vs)
{
    //  Get the corresponding ID used in switch statement
    unsigned long id = getLabel(vs[0]);
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    switch (id)
    {
            //  Construct global data
        case 0:     g_near    = toFloat(vs[1]);                 break;
        case 1:     g_left    = toFloat(vs[1]);                 break;
        case 2:     g_right   = toFloat(vs[1]);                 break;
        case 3:     g_bottom  = toFloat(vs[1]);                 break;
        case 4:     g_top     = toFloat(vs[1]);                 break;
        case 7:     g_back    = toVec4(vs[1], vs[2], vs[3]);    break;
        case 8:     g_ambient = toVec4(vs[1], vs[2], vs[3]);    break;
        case 9:     filename  = vs[1];                          break;
        case 5:
        {
            //  Construct sphere properties for each sphere
            Sphere sphere;
            sphere.name       = vs[1];
            sphere.pos        = toVec4(vs[2], vs[3], vs[4]);
            sphere.color      = toVec4(vs[8], vs[9], vs[10]);
            sphere.ambient    = toFloat(vs[11]);
            sphere.diffuse    = toFloat(vs[12]);
            sphere.specular   = toFloat(vs[13]);
            sphere.reflection = toFloat(vs[14]);
            sphere.specExp    = toFloat(vs[15]);
            sphere.trans      = Translate(sphere.pos) * Scale(toFloat(vs[5]),
                                                              toFloat(vs[6]),
                                                              toFloat(vs[7]));
            InvertMatrix(sphere.trans, sphere.inverseTrans);
            g_sphere.push_back(sphere);
            break;
        }
        case 6:
        {
            //  Construct light properties for each light source
            Light light;
            light.name      = vs[1];
            light.pos       = toVec4(vs[2], vs[3], vs[4]);
            light.intensity = toVec4(vs[5], vs[6], vs[7]);
            g_light.push_back(light);
            break;
        }
        case 10:    //  Get output size
            g_width = stoi(vs[1]);
            g_height = stoi(vs[2]);
            g_colors.resize(g_width * g_height);
            break;
        default:
            break;
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}

// -------------------------------------------------------------------
// Intersection routine

//  Return the calculated solutions to the quadratic function

void findSoln(float A, float B, float C, float delta, float& tMin, float& tMax)
{
    float t1 = -B / A;
    float variance = sqrt(delta) / A;
    float t2 = t1 + variance;
    t1 -= variance;

    tMin = min(t1, t2);
    tMax = max(t1, t2);
}

//  Get the normal for the transformed object

inline
vec4 getNormal(const vec3& start, const vec3& dir, float t, const Sphere& s)
{
    return normalize(vec4(toVec3(transpose(s.inverseTrans) * vec4(start + dir * t, 1)), 0));
}

//  Find out whether the specific ray intersects with a sphere
//  that satisfies requirements

bool intersectSphere(const vec4& start, const vec4& dir, vec4& normal,
                     float& t, const Sphere& s, bool traceLight, bool restriction)
{
    //  Inverse transform the starting point and the direction vector
    vec3 startP(toVec3(s.inverseTrans * start));
    vec3 dirP(toVec3(s.inverseTrans * dir));
    
    //  Calculate parameters for the quadratic equation
    float A = dot(dirP, dirP);
    float B = dot(startP, dirP);
    float C = dot(startP, startP) - 1;
    float delta = pow(B, 2) - A * C;
    
    //  Calculate the restriction
    const float minDist = (traceLight || ! restriction) ? RECOLLISION : MIN_DIST;
    
    //  One solution
    if (delta == 0.0f)
    {
        t = -B / A;
        if ( ! traceLight)
            normal = getNormal(startP, dirP, t, s);
        return t >= minDist;
    }
    
    //  No solution
    if (delta < 0.0f)
        return false;

    //  Two solutions
    float t1;
    findSoln(A, B, C, delta, t, t1);
    
    //  If the smaller t does not satisfy the requirement,
    //  Set it to the larger value.
    if (t < minDist)
        t = t1;
    if ( ! traceLight)
        normal = t == t1 ? -getNormal(startP, dirP, t, s) : getNormal(startP, dirP, t, s);
    return t >= minDist;    //  Check the requirement
}

//  Return the closest intersection from the input ray

bool intersect(const vec4& start, const vec4& dir, vec4& normal,
               float& t, Sphere& s, bool traceLight, bool restrcn)
{
    float tMin = -1, tPara;
    bool intersected = false;
    vec4 testNormal;
    
    //  Loop through all spheres
    for (auto it = g_sphere.begin(); it != g_sphere.end(); it++)
    {
        if (! intersectSphere(start, dir, testNormal, tPara, *it, traceLight, restrcn))
            continue;   //  No intersection, continue to the next available sphere
        //  Get the closest interesction
        if ( ! intersected || tPara < tMin)
        {
            normal = testNormal;
            tMin = tPara;
            s = *it;
            intersected = true;
        }
    }
    t = tMin;
    return intersected;
}

//  Wrapper for finding out shadow rays

inline
bool shadowRay(const vec4& start, const vec4& dir, float& t)
{
    vec4 dummy;
    Sphere s;
    bool intsec = intersect(start, dir, dummy, t, s, true, true);
    return intsec && t < 1;
}

//  Wrapper for finding out intersected spheres

inline
bool hasObject(const vec4& start, const vec4& dir, vec4& normal,
               float& t, Sphere& s, bool restrcn)
{
    return intersect(start, dir, normal, t, s, false, restrcn);
}

// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int depth)
{
    //  Color vectors with no contribution
    vec4 color(0.0f, 0.0f, 0.0f, 1.0f);
    
    //  Check Recursive Depth
    if (depth == 3)
        return color;
    
    //  Initialize local parameters
    vec4 start(ray.origin);
    vec4 dir(ray.dir);
    float t;
    Sphere s;
    vec4 normal;
    bool firstLevel = ! depth;

    //  Find closest intersection between objects
    if ( ! hasObject(start, dir, normal, t, s, firstLevel))
        return firstLevel ? g_back : color;

    //  Local Illumination

    //  Ambient
    color += s.ambient * g_ambient * s.color;

    //  Light Source
    const vec4 pointOfIntersect = vec4(toVec3(start + dir * t), 1);
    vec4 vision = -normalize(dir);

    for (auto it = g_light.begin(); it != g_light.end(); it++)
    {
        vec4 end = it->pos;
        //  (1-t) start + t (end) = start + t (end - start)
        vec4 light = end - pointOfIntersect;
        if (shadowRay(pointOfIntersect, light, t))
            continue;
        
        //  Normalize light direction pointers
        light = normalize(light);
        float dotProduct = dot(normal, light);
        
        //  Reflection vector
        vec4 reflection = normalize(2 * dotProduct * normal - light);
        
        //  Local Illumination Formula
        color += (it->intensity * max(dotProduct, 0.0f) * s.color * s.diffuse
                + it->intensity * powf(max(dot(reflection, vision), 0.0f), s.specExp) * s.specular);
    }

    //  Reflection

    Ray ref;
    ref.origin = pointOfIntersect;
    ref.dir = 2 * dot(normal, vision) * normal - vision;
    
    //  Findout reflections recursively
    color += s.reflection * trace(ref, depth + 1);

    //  Return results
    return normalizeColorComponent(color);
}

//  Calculate resulting point in getDir function
inline
float calcParameter(float start, float end, float alpha)
{
    return alpha * end + (1.0 - alpha) * start;
}

//  Get directions for ray tracing

vec4 getDir(int ix, int iy)
{
    //  Do not normalize at this point to
    //  get easier checkings for sphere cullings
    return vec4(calcParameter(g_left, g_right, (float) ix / g_width),
                calcParameter(g_bottom, g_top, (float) iy / g_height),
                -g_near, 0.0f);
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    int depth = 0;
    vec4 color = trace(ray, depth);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;
    
    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);
    
    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }
    
    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
        {
            int colorNum = y * g_width + x;
            int buckNum = 3 * colorNum;
            buf[buckNum]   = (unsigned char) (g_colors[colorNum].x * 255.9f);
            buf[buckNum+1] = (unsigned char) (g_colors[colorNum].y * 255.9f);
            buf[buckNum+2] = (unsigned char) (g_colors[colorNum].z * 255.9f);
        }
    
    //  Output with the correct filename
    savePPM(g_width, g_height, filename.c_str(), buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: raytrace <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
    return 0;
}

