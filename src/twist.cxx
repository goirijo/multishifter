#include "./twist.hpp"
#include "./misc.hpp"
#include <filesystem>
#include <casmutils/mush/twist.hpp>
#include <multishift/slice_settings.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <ostream>
#include <utility>
#include <vector>

std::string make_twist_dirname(double cleavage)
{
    std::stringstream twiststream;
    twiststream << std::fixed << std::setprecision(6) << cleavage;
    return "twist__" + twiststream.str();
}

void setup_subcommand_twist(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* twist_sub = app.add_subcommand("twist", "Create approximated supercells that can accommodate moirons on a rotated bilayer.");
    twist_sub->add_option("-s,--settings", *settings_path_ptr, "Settings file with slab, and list of angles to rotate along.");

    twist_sub->callback([settings_path_ptr]() { run_subcommand_twist(*settings_path_ptr, std::cout); });
}

/* void write_twisted_structures(const std::vector<mush::MoireStructureGenerator>& generators, const mush::fs::path& root, std::ostream& log) */
/* { */
/*     using LATTICE=mush::MoireStructureGenerator::LATTICE; */
/*     std::pair<LATTICE,std::string> A{LATTICE::ALIGNED,"aligned"}; */
/*     std::pair<LATTICE,std::string> R{LATTICE::ROTATED,"rotated"}; */

/*     log<<"Write cleaved structures to..."<<root<<std::endl; */
/*     for(const auto& gen : generators) */
/*     { */
/*         auto target=root/make_twist_dirname(gen.degrees()); */
/*         cautious_create_directory(target); */

/*         for(auto zone_pr : {A,R}) */
/*         { */
/*             for(auto lat_pr : {A,R}) */
/*             { */
/*                 auto layer=gen.layer(zone_pr.first,lat_pr.first); */
/*                 std::string combo_name=zone_pr.second+"_zone_"+lat_pr.second+"_lattice.vasp"; */
/*                 mush::cu::xtal::write_poscar(layer,target/combo_name); */
/*             } */
/*         } */
/*     } */
/* } */

void run_subcommand_twist(const mush::fs::path& settings_path, std::ostream& log)
{
    /* json settings = load_json(settings_path); */
    /* std::string project_name = extract_name_from_settings(settings); */

    /* auto slaber_settings = mush::SlabSettings::from_json(settings); */
    /* auto twist_settings = mush::TwisterSettings::from_json(settings); */

    /* log << "Reading aligned surface from " << slaber_settings.slab_unit_path << "...\n"; */
    /* auto slab_unit = mush::cu::xtal::Structure::from_poscar(slaber_settings.slab_unit_path); */

    /* log << "Twisting layer..."; */
    /* std::vector<mush::MoireStructureGenerator> all_twists; */
    /* for(double angle : twist_settings.angles) */
    /* { */
    /*     log<<".."; */
    /*     all_twists.emplace_back(slab_unit,angle); */
    /* } */
    /* log << "\n"; */

    /* mush::fs::path twisted_path(project_name+".twist"); */
    /* cautious_create_directory(twisted_path); */
    /* write_twisted_structures(all_twists,twisted_path,log); */

    /* log << "Back up settings used...\n"; */
    /* write_json(settings, twisted_path / "twist.json"); */

    /* log << "Back up slab...\n"; */
    /* mush::cu::xtal::write_poscar(slab_unit,twisted_path/ "slab.vasp"); */


    return;
}
