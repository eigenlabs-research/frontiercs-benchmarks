#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace std::chrono;

struct Point3D { double x, y, z; };

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int n;
    if (!(cin >> n)) return 0;
    
    auto start_time = high_resolution_clock::now();
    
    int m = ceil(cbrt(n));
    vector<Point3D> pts;
    pts.reserve(n);
    
    double step = (m > 1) ? 1.0 / (m - 1) : 0.5;
    double offset = (m > 1) ? 0.0 : 0.5;
    
    for(int i=0; i<m && pts.size() < n; ++i) {
        for(int j=0; j<m && pts.size() < n; ++j) {
            for(int k=0; k<m && pts.size() < n; ++k) {
                pts.push_back({i*step + offset, j*step + offset, k*step + offset});
            }
        }
    }
    
    double margin = 1.0 / (2.0 * m);
    for(auto& p : pts) {
        p.x = margin + p.x * (1.0 - 2.0*margin);
        p.y = margin + p.y * (1.0 - 2.0*margin);
        p.z = margin + p.z * (1.0 - 2.0*margin);
    }
    
    double step_size = 0.005;
    int iterations = 0;
    
    while(true) {
        if ((iterations & 7) == 0) {
            auto now = high_resolution_clock::now();
            if (duration_cast<milliseconds>(now - start_time).count() > 1850) break;
        }
        iterations++;
        
        vector<Point3D> forces(n, {0,0,0});
        double min_d2 = 1e9;
        
        for(int i=0; i<n; ++i) {
            for(int j=i+1; j<n; ++j) {
                double dx = pts[i].x - pts[j].x;
                double dy = pts[i].y - pts[j].y;
                double dz = pts[i].z - pts[j].z;
                double d2 = dx*dx + dy*dy + dz*dz;
                
                if (d2 < min_d2) min_d2 = d2;
                
                if (d2 < 0.1 && d2 > 1e-9) {
                    double f = 1.0 / (d2 * sqrt(d2));
                    forces[i].x += f * dx; forces[i].y += f * dy; forces[i].z += f * dz;
                    forces[j].x -= f * dx; forces[j].y -= f * dy; forces[j].z -= f * dz;
                }
            }
            
            double d_x0 = pts[i].x; if (d_x0 < 0.2) forces[i].x += 1.0 / (d_x0*d_x0);
            double d_x1 = 1.0 - pts[i].x; if (d_x1 < 0.2) forces[i].x -= 1.0 / (d_x1*d_x1);
            
            double d_y0 = pts[i].y; if (d_y0 < 0.2) forces[i].y += 1.0 / (d_y0*d_y0);
            double d_y1 = 1.0 - pts[i].y; if (d_y1 < 0.2) forces[i].y -= 1.0 / (d_y1*d_y1);
            
            double d_z0 = pts[i].z; if (d_z0 < 0.2) forces[i].z += 1.0 / (d_z0*d_z0);
            double d_z1 = 1.0 - pts[i].z; if (d_z1 < 0.2) forces[i].z -= 1.0 / (d_z1*d_z1);
        }
        
        double max_f = 0;
        for(int i=0; i<n; ++i) {
            double fm = sqrt(forces[i].x*forces[i].x + forces[i].y*forces[i].y + forces[i].z*forces[i].z);
            if (fm > max_f) max_f = fm;
        }
        
        if (max_f > 1e-9) {
            for(int i=0; i<n; ++i) {
                pts[i].x += step_size * (forces[i].x / max_f);
                pts[i].y += step_size * (forces[i].y / max_f);
                pts[i].z += step_size * (forces[i].z / max_f);
                
                if (pts[i].x < 0.0) pts[i].x = 0.0;
                if (pts[i].x > 1.0) pts[i].x = 1.0;
                if (pts[i].y < 0.0) pts[i].y = 0.0;
                if (pts[i].y > 1.0) pts[i].y = 1.0;
                if (pts[i].z < 0.0) pts[i].z = 0.0;
                if (pts[i].z > 1.0) pts[i].z = 1.0;
            }
        }
        
        step_size *= 0.999;
    }
    
    for(int i=0; i<n; ++i) {
        cout << pts[i].x << " " << pts[i].y << " " << pts[i].z << "\n";
    }
    
    return 0;
}
