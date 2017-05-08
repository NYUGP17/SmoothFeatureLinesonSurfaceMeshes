//
// Created by Liang Zhuo on 5/8/17.
//

#ifndef SMOOTHFEATURELINE_UFSET_H
#define SMOOTHFEATURELINE_UFSET_H


class UFSet {
public:
    UFSet(int);
    int find(int);
    void merge(int, int);
protected:
    std::vector <int> s;
};


#endif //SMOOTHFEATURELINE_UFSET_H
