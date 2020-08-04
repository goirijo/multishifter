#include "./translate.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "misc.hpp"
#include "multishift/slab.hpp"
#include <filesystem>
#include <multishift/shift.hpp>
#include <multishift/slice_settings.hpp>
#include <multishift/twist.hpp> //defines frankenstein

std::string make_translate_dirname(int a, int b) { return "translate__" + std::to_string(a) + "." + std::to_string(b); }

void setup_subcommand_translate(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* translate_sub = app.add_subcommand("translate", "Translate the basis of a structure by in plane vectors of a grid.");
    translate_sub
        ->add_option("-s,--settings", *settings_path_ptr, "Settings file with path to slab unit, slab thickness, and grid density.")
        ->required();

    translate_sub->callback([settings_path_ptr]() { run_subcommand_translate(*settings_path_ptr, std::cout); });
}

void run_subcommand_translate(const mush::fs::path& settings_path, std::ostream& log)
{
    json settings = load_json(settings_path);
    std::string project_name = extract_name_from_settings(settings);
    mush::fs::path translate_root(project_name + ".translate");

    log << "Project name: " << project_name << std::endl;
    cautious_create_directory(translate_root);

    auto translate_settings = mush::ShiftSettings::from_json(settings);
    auto slaber_settings = mush::SlabSettings::from_json(settings);
    auto slab = mush::make_stacked_slab(Structure::from_poscar(slaber_settings.slab_unit_path), slaber_settings.stacks);

    auto [translations, records] =
        mush::make_uniform_in_plane_shift_vectors(slab.lattice(), translate_settings.a_density, translate_settings.b_density);

    assert(translations.size() == records.size());

    /* json final_record; */
    log << "Translating structures...\n";
    for (int i = 0; i < translations.size(); ++i)
    {

        const auto& sr = records[i];
        /* json sr_json; */
        /* sr_json["a_index"] = sr.a; */
        /* sr_json["b_index"] = sr.b; */

        mush::fs::path ttarget = translate_root / make_translate_dirname(sr.a, sr.b);
        /* final_record[ttarget]=sr_json; */

        auto translated_struc = mush::cu::xtal::frankenstein::translate_basis(slab, translations[i]);

        cautious_create_directory(ttarget);
        mush::cu::xtal::write_poscar(translated_struc, ttarget / "POSCAR");
    }

    log << "Wrote translated structures to..." << translate_root << std::endl;

    /* log << "Save record of structures to " << translate_root << "\n"; */
    /* write_json(final_record,translate_root/"record.json"); */

    log << "Save record of slab " << translate_root << "\n";
    mush::cu::xtal::write_poscar(slab,translate_root/"slab.vasp");

    log << "Back up settings used...\n";
    write_json(settings,translate_root/"settings.json");

    return;
}
