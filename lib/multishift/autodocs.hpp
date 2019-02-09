#ifndef MULTIDOCS_HH
#define MULTIDOCS_HH

#include <string>
#include <tuple>
#include <unordered_map>
#include "casm/casm_io/jsonParser.hh"

namespace mush
{
    namespace docs
    {
        /**
         * Enum to sanely index into the tuples that the settings info
         * are stored as
         */

        enum class SettingSlotIndex
        {
            ENTRY,
            DESC,
            REQUIRED,
            DEFAULT
        };

        /**
         * This class is meant to keep documentation in a single place,
         * and provide an interface with which to extract desctiptions
         * for different settings, and automatically construct example
         * json settings.
         */

        class SettingsInfo
        {
            public:

                ///Collection of information for a particular setting
                typedef std::tuple<std::string,std::string,bool,std::string> SettingSlot;
                typedef SettingSlotIndex SSI;

                ///Initialize with the tag that the entire settings are under within the json
                SettingsInfo(const std::string& init_tag);

                const std::string& tag() const;

                ///Program options style chaining of settings options, short description, whether it's required,
                ///and default value.
                SettingsInfo& add_setting();

                ///Documentation string for each of the entries of the setting
                std::unordered_map<std::string,std::string> descriptions() const;
                
                ///Example json format
                CASM::jsonParser example() const;

            private:

                std::string m_tag;
        };
    }
}

#endif
