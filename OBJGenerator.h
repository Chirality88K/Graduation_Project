#ifndef OBJGENERATOR
#define OBJGENERATOR
#include <vector>
#include <string>
#include <fstream>
#include "write3dm.h"
#include "ParameterSurface.h"

class OBJData
{
public:
    OBJData() {}
    OBJData(const OBJData &data) : name_(data.name_), points_(data.points_), normals_(data.normals_), faces_(data.faces_) {}
    void SetName(const std::string &name)
    {
        name_ = name;
    }
    std::string GetName() const
    {
        return name_;
    }
    unsigned int GetVertexNumber() const
    {
        return points_.size() / 3;
    }
    unsigned int AddPoint(double px, double py, double pz);
    unsigned int AddNormal(double nx, double ny, double nz);
    unsigned int AddQuadFace(unsigned int index1, unsigned int index2, unsigned int index3, unsigned int index4);
    unsigned int AddTriangleFace(unsigned int index1, unsigned int index2, unsigned int index3);
    void WritePoints(std::ofstream &filestream) const;
    void WriteNormals(std::ofstream& filestream) const;
    void WriteFace_with_IndexOffset(std::ofstream& filestream, unsigned int offset = 0) const;

private:
    std::string name_ = "Default";
    std::vector<double> points_;
    std::vector<double> normals_;
    std::vector<unsigned int> faces_;
};

class OBJGenerator
{
public:
    OBJGenerator() {}
    ~OBJGenerator()
    {
        for (OBJData *pdata : objdata_)
        {
            delete pdata;
        }
    }
    void AddData_Copy(const OBJData &data);
    void AddData_NoCopy(OBJData *data);
    void AddNurbsSurface(const ON_NurbsSurface& ons, const std::string& name, unsigned int u_sample_knot, unsigned int v_sample_count, bool u_period = false, bool v_period = false);
    void AddParameterSurface(const ParameterSurface& ps, const std::string& name, unsigned int u_sample_knot, unsigned int v_sample_count, bool u_period = false, bool v_period = false);
    bool Write(const std::string &filename) const;

private:
    std::vector<OBJData *> objdata_;
};

#endif