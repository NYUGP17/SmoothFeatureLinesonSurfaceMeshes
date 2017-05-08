//
// Created by Liang Zhuo on 5/8/17.
//

#include <vector>

#include "UFSet.h"

UFSet::UFSet(int n) {
    s = std::vector <int>(n, -1);
}

int UFSet::find(int x) {
    return (s[x] < 0) ? x : (s[x] = find(s[x]));
}

void UFSet::merge(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y)
        return;
    if (s[x] > s[y])
        std::swap(x, y);
    else if (s[x] == s[y])
        -- s[x];
    s[y] = x;
}