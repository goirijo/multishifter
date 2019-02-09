#include "./autodocs.hpp"

namespace mush
{
    namespace docs
    {
        SettingsInfo::SettingsInfo(const std::string& init_tag):m_tag(init_tag){}
        const std::string& SettingsInfo::tag() const {return m_tag;};
    }
}
