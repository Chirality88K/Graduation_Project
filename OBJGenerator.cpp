#include "OBJGenerator.h"

extern std::string output_base_dir;

unsigned int OBJData::AddPoint(double px, double py, double pz)
{
    unsigned int N = points_.size() / 3;
    points_.push_back(px);
    points_.push_back(py);
    points_.push_back(pz);
    return N + 1;
}

unsigned int OBJData::AddNormal(double nx, double ny, double nz)
{
    unsigned int N = normals_.size() / 3;
    double l = sqrt(pow(nx, 2) + pow(ny, 2) + pow(nz, 2));
    normals_.push_back(nx / l);
    normals_.push_back(ny / l);
    normals_.push_back(nz / l);
    return N + 1;
}

unsigned int OBJData::AddQuadFace(unsigned int index1, unsigned int index2, unsigned int index3, unsigned int index4)
{
    AddTriangleFace(index1, index2, index3);
    return AddTriangleFace(index3, index4, index1);
}

unsigned int OBJData::AddTriangleFace(unsigned int index1, unsigned int index2, unsigned int index3)
{
    unsigned int N = faces_.size() / 3;
    faces_.push_back(index1);
    faces_.push_back(index2);
    faces_.push_back(index3);
    return N + 1;
}

void OBJData::WritePoints(std::ofstream &fs) const
{
    fs << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < points_.size(); i = i + 3)
    {
        fs << "v " << points_[i] << " " << points_[i + 1] << " " << points_[i + 2] << "\n";
    }
}

void OBJData::WriteNormals(std::ofstream& fs) const
{
    fs << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < normals_.size(); i = i + 3)
    {
        fs << "vn " << normals_[i] << " " << normals_[i + 1] << " " << normals_[i + 2] << "\n";
    }
}

void OBJData::WriteFace_with_IndexOffset(std::ofstream& fs, unsigned int offset) const
{
    for (size_t i = 0; i < faces_.size(); i = i + 3)
    {
        fs << "f " <<
            faces_[i + 0] + offset << "/1/" << faces_[i + 0] + offset << " " <<
            faces_[i + 1] + offset << "/1/" << faces_[i + 1] + offset << " " <<
            faces_[i + 2] + offset << "/1/" << faces_[i + 2] + offset << " " <<
            "\n";
    }
}

void OBJGenerator::AddData_Copy(const OBJData &data)
{
    objdata_.push_back(new OBJData(data));
}

void OBJGenerator::AddData_NoCopy(OBJData *data)
{
    objdata_.push_back(data);
}

void OBJGenerator::AddNurbsSurface(const ON_NurbsSurface& ons, const std::string& name, unsigned int u_sample_count, unsigned int v_sample_count, bool u_period, bool v_period)
{
    OBJData* data = new OBJData();
    data->SetName(name);
    double u0, u1, v0, v1;
    ons.GetDomain(0, &u0, &u1);
    ons.GetDomain(1, &v0, &v1);
    unsigned int u_count = u_period ? u_sample_count : u_sample_count + 1;
    unsigned int v_count = v_period ? v_sample_count : v_sample_count + 1;
    for (int i = 0; i < u_count; ++i)
    {
        for (int j = 0; j < v_count; ++j)
        {
            double u = (u1 - u0) / u_sample_count * i + u0;
            double v = (v1 - v0) / v_sample_count * j + v0;
            ON_3dPoint p = ons.PointAt(u, v);
            ON_3dVector n = ons.NormalAt(u, v);
            data->AddPoint(p.x, p.y, p.z);
            data->AddNormal(n.x, n.y, n.z);
        }
    }
    unsigned int u_f_cnt = u_sample_count;
    unsigned int v_f_cnt = v_sample_count;
    auto point_index = [u_count,v_count](int i, int j)->int {
        i = i % u_count;
        j = j % v_count;
        return i * v_count + j + 1;
    };

    for (int i = 0; i < u_f_cnt; ++i)
    {
        for (int j = 0; j < v_f_cnt; ++j)
        {
            data->AddQuadFace(point_index(i, j), point_index(i + 1, j), point_index(i + 1, j + 1), point_index(i, j + 1));
        }
    }
    this->AddData_NoCopy(data);
}

void OBJGenerator::AddParameterSurface(const ParameterSurface& ps, const std::string& name, unsigned int u_sample_count, unsigned int v_sample_count, bool u_period, bool v_period)
{
    OBJData* data = new OBJData();
    data->SetName(name);
    double u0, u1, v0, v1;
    ps.GetDomain(0, &u0, &u1);
    ps.GetDomain(1, &v0, &v1);
    unsigned int u_count = u_period ? u_sample_count : u_sample_count + 1;
    unsigned int v_count = v_period ? v_sample_count : v_sample_count + 1;
    for (int i = 0; i < u_count; ++i)
    {
        for (int j = 0; j < v_count; ++j)
        {
            double u = (u1 - u0) / u_sample_count * i + u0;
            double v = (v1 - v0) / v_sample_count * j + v0;
            ON_3dPoint p = ps.PointAt(u, v);
            ON_3dVector n = ps.NormalWithoutCheckUnit(u, v);
            data->AddPoint(p.x, p.y, p.z);
            data->AddNormal(n.x, n.y, n.z);
        }
    }
    unsigned int u_f_cnt = u_sample_count;
    unsigned int v_f_cnt = v_sample_count;
    auto point_index = [u_count, v_count](int i, int j)->int {
        i = i % u_count;
        j = j % v_count;
        return i * v_count + j + 1;
    };

    for (int i = 0; i < u_f_cnt; ++i)
    {
        for (int j = 0; j < v_f_cnt; ++j)
        {
            data->AddQuadFace(point_index(i, j), point_index(i + 1, j), point_index(i + 1, j + 1), point_index(i, j + 1));
        }
    }
    this->AddData_NoCopy(data);
}

bool OBJGenerator::Write(const std::string &filename_without_extension) const
{
    std::string filename = output_base_dir + filename_without_extension + ".obj";
    std::ofstream ofs(filename);
    unsigned int index_offset = 0;
    for (OBJData *pdata : objdata_)
    {
        ofs << "o " << pdata->GetName() << "\n";
        pdata->WritePoints(ofs);
        ofs << "vt 0 0\n";
        pdata->WriteNormals(ofs);
        pdata->WriteFace_with_IndexOffset(ofs, index_offset);
        index_offset += pdata->GetVertexNumber();
    }
    ofs.close();
    return true;
}
