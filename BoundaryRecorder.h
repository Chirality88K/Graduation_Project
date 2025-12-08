#ifndef BOUNDARYRECORDER_H
#define BOUNDARYRECORDER_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <string>
#include <vector>

class Boundary
{
public:
    Boundary(const ON_3dPoint &ps, const ON_3dPoint &pe, const ON_3dVector &vs, const ON_3dVector &ve) : ps_(ps), pe_(pe), vs_(vs), ve_(ve) {}

public:
    ON_3dPoint ps_;
    ON_3dPoint pe_;
    ON_3dVector vs_;
    ON_3dVector ve_;
};

class BoundaryRecorder
{
public:
    static BoundaryRecorder& GetRecorder()
    {
        static BoundaryRecorder br_;
        return br_;
    }
    void Append(const ON_3dPoint &ps, const ON_3dPoint &pe, const ON_3dVector &vs, const ON_3dVector &ve)
    {
        boundary_.push_back(Boundary(ps, pe, vs, ve));
    }
    void Write() const;
    void Read();
    const std::vector<Boundary>& GetBoundary() const
    {
        return boundary_;
    }

private:
    std::vector<Boundary> boundary_;
private:
    BoundaryRecorder() {}
    BoundaryRecorder(const BoundaryRecorder&) = delete;
    BoundaryRecorder(BoundaryRecorder&&) = delete;
    BoundaryRecorder& operator=(const BoundaryRecorder&) = delete;
    BoundaryRecorder& operator=(BoundaryRecorder&&) = delete;
};

#endif