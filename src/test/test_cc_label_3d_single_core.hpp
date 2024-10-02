//
// Created by muyezhu on 12/7/17.
//

#ifndef MCP3D_TEST_CC_LABEL_3D_SINGLE_CORE_HPP
#define MCP3D_TEST_CC_LABEL_3D_SINGLE_CORE_HPP
namespace mcp3d
{
namespace test
{
void NdiLabel();

bool LabelCC3DCorrect();

bool TestEqLabelSet();

template <typename T = int32>
bool CompareWithNDILabel(T* label, int32* ndi_label)
{
    return false;
}

}
}

#endif //MCP3D_TEST_CC_LABEL_3D_SINGLE_CORE_HPP
