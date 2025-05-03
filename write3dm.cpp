#include "write3dm.h"

void ChiralityWrite3dmModel(const ONX_Model* model, const wchar_t* filename)
{
    ON_TextLog error_log;
    bool success = model->Write(filename, 0, &error_log);
    if (success)
    {
        std::cout << "Successfully wrote ";
        std::wcout << std::wstring(filename) << "\n";
    }
    else
    {
        std::cout << "Fail to wrtie ";
        std::wcout << std::wstring(filename) << "\n";
    }
}

ON_3dmObjectAttributes* Internal_CreateManagedAttributes(int layer_index, const wchar_t* name)
{
    ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
    attributes->m_layer_index = layer_index;
    attributes->m_name = name;
    return attributes;
}