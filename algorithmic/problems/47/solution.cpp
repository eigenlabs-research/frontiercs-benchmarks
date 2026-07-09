#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace std::chrono;

struct Item {
    string type;
    int w, h, v, limit;
};

struct Placement {
    string type;
    int x, y, rot;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    string input_data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    if (input_data.empty()) return 0;
    
    int W = 0, H = 0;
    bool allow_rotate = false;
    
    smatch m;
    if (regex_search(input_data, m, regex("\"W\"\\s*:\\s*(\\d+)"))) W = stoi(m[1].str());
    if (regex_search(input_data, m, regex("\"H\"\\s*:\\s*(\\d+)"))) H = stoi(m[1].str());
    if (regex_search(input_data, m, regex("\"allow_rotate\"\\s*:\\s*(true|false)"))) allow_rotate = (m[1].str() == "true");
    
    vector<Item> items;
    regex item_regex("\\{[^\\}]*\"type\"\\s*:\\s*\"([^\"]+)\"[^\\}]*\\}");
    
    auto items_begin = sregex_iterator(input_data.begin(), input_data.end(), item_regex);
    auto items_end = sregex_iterator();
    
    for (sregex_iterator i = items_begin; i != items_end; ++i) {
        string item_str = i->str();
        Item it;
        smatch m2;
        if (regex_search(item_str, m2, regex("\"type\"\\s*:\\s*\"([^\"]+)\""))) it.type = m2[1].str();
        if (regex_search(item_str, m2, regex("\"w\"\\s*:\\s*(\\d+)"))) it.w = stoi(m2[1].str());
        if (regex_search(item_str, m2, regex("\"h\"\\s*:\\s*(\\d+)"))) it.h = stoi(m2[1].str());
        if (regex_search(item_str, m2, regex("\"v\"\\s*:\\s*(\\d+)"))) it.v = stoi(m2[1].str());
        if (regex_search(item_str, m2, regex("\"limit\"\\s*:\\s*(\\d+)"))) it.limit = stoi(m2[1].str());
        items.push_back(it);
    }
    
    sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        double d1 = (double)a.v / (a.w * a.h);
        double d2 = (double)b.v / (b.w * b.h);
        return d1 > d2;
    });
    
    vector<Placement> placements;
    vector<int> counts(items.size(), 0);
    
    int current_y = 0;
    int current_x = 0;
    int shelf_h = 0;
    
    while(true) {
        bool placed_any = false;
        for (size_t i = 0; i < items.size(); ++i) {
            if (counts[i] >= items[i].limit) continue;
            
            int pw = items[i].w;
            int ph = items[i].h;
            
            if (current_x + pw <= W && current_y + ph <= H) {
                placements.push_back({items[i].type, current_x, current_y, 0});
                current_x += pw;
                shelf_h = max(shelf_h, ph);
                counts[i]++;
                placed_any = true;
                break;
            }
            
            if (allow_rotate) {
                pw = items[i].h;
                ph = items[i].w;
                if (current_x + pw <= W && current_y + ph <= H) {
                    placements.push_back({items[i].type, current_x, current_y, 1});
                    current_x += pw;
                    shelf_h = max(shelf_h, ph);
                    counts[i]++;
                    placed_any = true;
                    break;
                }
            }
        }
        
        if (!placed_any) {
            if (shelf_h == 0) break; 
            current_y += shelf_h;
            current_x = 0;
            shelf_h = 0;
            if (current_y >= H) break;
            
            bool can_fit = false;
            for (size_t i = 0; i < items.size(); ++i) {
                if (counts[i] < items[i].limit) {
                    if (items[i].w <= W && current_y + items[i].h <= H) can_fit = true;
                    if (allow_rotate && items[i].h <= W && current_y + items[i].w <= H) can_fit = true;
                }
            }
            if (!can_fit) break;
        }
    }
    
    cout << "{\n  \"placements\": [\n";
    for(size_t i=0; i<placements.size(); ++i) {
        cout << "    {\"type\":\"" << placements[i].type << "\", \"x\":" << placements[i].x << ", \"y\":" << placements[i].y << ", \"rot\":" << placements[i].rot << "}";
        if (i + 1 < placements.size()) cout << ",";
        cout << "\n";
    }
    cout << "  ]\n}\n";
    
    return 0;
}
