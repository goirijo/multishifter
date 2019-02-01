#ifndef MULTISHIFTSETTINGS_HH
#define MULTISHIFTSETTINGS_HH

#include "./base.hpp"
#include "./exceptions.hpp"
#include "./shift.hpp"
#include <memory>

namespace mush
{

class FullSettings
{
public:
    static FullSettings from_json(const CASM::jsonParser& init_json);
    static FullSettings from_path(const fs::path& init_path);

    const std::string& name() const { return m_name; }

    const BaseSettings& base_settings() const;
    const ShiftSettings& shift_settings() const;

private:
    FullSettings(const std::string& init_name, std::unique_ptr<BaseSettings>&& init_base_settings,
                 std::unique_ptr<ShiftSettings>&& init_shift_settings);

    /// Name for the particular set of shifts you're going to do
    std::string m_name;

    /// Specifies all the settings involved with creating a slab
    std::unique_ptr<BaseSettings> m_base_settings_ptr;

    /// Specifies the settings for the gamma surface density and cleavage between slabs
    std::unique_ptr<ShiftSettings> m_shift_settings_ptr;
};

FullSettings::FullSettings(const std::string& init_name, std::unique_ptr<BaseSettings>&& init_base_settings,
                           std::unique_ptr<ShiftSettings>&& init_shift_settings)
    : m_name(init_name),
      m_base_settings_ptr(std::move(init_base_settings)),
      m_shift_settings_ptr(std::move(init_shift_settings))
{
}

FullSettings FullSettings::from_json(const CASM::jsonParser& init_json)
{
    std::unique_ptr<BaseSettings> base_settings_ptr(nullptr);
    std::unique_ptr<ShiftSettings> shift_settings_ptr(nullptr);

    if(!init_json.contains("name"))
    {
        throw except::UnspecifiedSettings("name");
    }

    if(init_json.contains(BaseSettings::tag()))
    {
        base_settings_ptr=std::unique_ptr<BaseSettings>(new BaseSettings(BaseSettings::from_json(init_json[BaseSettings::tag()])));
    }

    if(init_json.contains(ShiftSettings::tag()))
    {
        shift_settings_ptr=std::unique_ptr<ShiftSettings>(new ShiftSettings(ShiftSettings::from_json(init_json[ShiftSettings::tag()])));
    }

    return FullSettings(init_json["name"].get<std::string>(), std::move(base_settings_ptr), std::move(shift_settings_ptr));
}

FullSettings FullSettings::from_path(const fs::path& init_path)
{
    CASM::jsonParser settings_dump(init_path);
    return FullSettings::from_json(settings_dump);
}

const BaseSettings& FullSettings::base_settings() const
{
    if (m_base_settings_ptr.get() != nullptr)
    {
        return *m_base_settings_ptr.get();
    }
    else
    {
        throw except::UnspecifiedSettings(BaseSettings::tag());
    }
}

const ShiftSettings& FullSettings::shift_settings() const
{
    if (m_shift_settings_ptr.get() != nullptr)
    {
        return *m_shift_settings_ptr.get();
    }
    else
    {
        throw except::UnspecifiedSettings(ShiftSettings::tag());
    }
}
} // namespace mush

#endif
