#include "./twist.hpp"
#include "./common_options.hpp"
#include "./misc.hpp"
#include <casmutils/mush/twist.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <filesystem>
#include <ostream>
#include <utility>
#include <vector>

namespace cu = casmutils;

void setup_subcommand_twist(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto angles_ptr = std::make_shared<std::vector<double>>();
    auto max_lattice_sites_ptr = std::make_shared<int>();
    auto error_tol_ptr = std::make_shared<double>();
    auto zone_ptr = std::make_shared<std::string>();
    auto supercells_ptr = std::make_shared<std::string>();

    CLI::App* twist_sub =
        app.add_subcommand("twist", "Create approximate supercells that can accommodate emerging moirons from a specified rotation angle.");

    populate_subcommand_input_option(twist_sub, input_path_ptr.get());
    populate_subcommand_output_option(twist_sub, output_path_ptr.get());

    // clang-format off
    twist_sub->add_option("-a,--angles", *angles_ptr, "Rotation angles to twist the structure with in degrees. Rotation is applied at the origin, perpendicular to the ab-plane.")->required();
    twist_sub->add_option("-m,--max-lattice-sites", *max_lattice_sites_ptr, "Sets the maximum search space for more commensurate Moire supercells. If zero, don't try looking for supercells.")->default_val(0);
    twist_sub->add_option("-e,--error-tol", *error_tol_ptr, "Minimum improvement necessary to consider a larger supercell better than a smaller one.")->default_val(1e-8);
    twist_sub->add_option("-z,--brillouin-zone", *zone_ptr, "Which Brillouin zone to use when mapping reciprocal Moire lattice vectors back into the first Brillouin zone.")->default_val("aligned")->check(CLI::IsMember({"aligned","rotated"},CLI::ignore_case));
    twist_sub->add_option("-s,--supercells", *supercells_ptr, "Specify whether only the best supercell or every possible supercell that can hold max-lattice-sites should be saved.")->default_val("best")->check(CLI::IsMember({"best","all"},CLI::ignore_case));
    // clang-format off

    twist_sub->callback([=]() {run_subcommand_twist(*input_path_ptr,
            *output_path_ptr,
            *angles_ptr,
            *max_lattice_sites_ptr,
            *error_tol_ptr,
            *zone_ptr,
            *supercells_ptr,
            std::cout); });
}

std::string lat_to_name(mush::MoireStructureReport::LATTICE lat)
{
    return lat==mush::MoireStructureReport::LATTICE::ALIGNED ? "top" : "bottom";
}

std::string zone_to_name(mush::MoireStructureReport::ZONE zone)
{
    return zone==mush::MoireStructureReport::ZONE::ALIGNED ? "aligned" : "rotated";
}

std::string make_twist_dirname(double twist)
{
    std::stringstream twiststream;
    twiststream << std::fixed << std::setprecision(6) << twist;
    return "twist__" + twiststream.str();
}

std::string make_twist_id(double twist, const mush::MoireLatticeReport& report)
{
    std::stringstream twiststream;
    twiststream << std::fixed << std::setprecision(6) << twist;

    twiststream<<":"<<report.num_moirons();
    return "t" + twiststream.str();
}

mush::fs::path make_target_structure_dir(double twist, mush::MoireStructureReport::LATTICE lat, int scel_ix)
{
    if(scel_ix<0)
    {
        return make_twist_dirname(twist);
    }

    return mush::fs::path(make_twist_dirname(twist))/std::to_string(scel_ix);
}


namespace Eigen
{
    using mush::json;

    template<typename MatType>
    void to_json(json& j, const MatType& mat)
    {
        std::vector<std::vector<double>> casted;
        for(int i=0; i<mat.rows(); ++i)
        {
            const auto r=mat.row(i);
            std::vector<double> c;
            for(int j=0; j<r.size(); ++j)
            {
                c.push_back(r(j));
            }
            casted.emplace_back(c);
        }

        j=mush::json(casted);
    }
}

mush::json serialize(const mush::MoireStructureReport& moire)
{
    mush::json report;

    report["brillouin_zone"]=zone_to_name(moire.zone);
    /* report["true_moire_supercell_matrix"]=moire.true_moire_supercell_matrix; */
    /* ... */

    mush::DeformationReport def_report(moire.approximation_deformation);
    report["deformation"]=def_report.deformation;
    report["rotation"]=def_report.rotation;
    report["strain"]=def_report.strain;
    report["angle_adjustment"]=def_report.rotation_angle;
    report["adjustment_strain_metrics"]=def_report.strain_metrics;
    report["dilation_adjustment"]=def_report.dilation_strain;
    report["deviatoric_adjustment"]=def_report.deviatoric_strain;

    return report;
}

void commit_twisted_id(double twist, const mush::MoireStructureReport& best_report, mush::json* _record, const mush::fs::path& output_dir, const mush::fs::path& root, mush::MoireStructureReport::LATTICE lat)
{
    mush::json& twist_record=*_record;

                auto id=make_twist_id(twist,best_report);
                
                twist_record[id][lat_to_name(lat)]=serialize(best_report);
                auto tile_path=root/(lat_to_name(lat)+"_tile.vasp");
                twist_record[id][lat_to_name(lat)]["tile"]=tile_path;
                auto layer_path=root/(lat_to_name(lat)+"_layer.vasp");
                twist_record[id][lat_to_name(lat)]["layer"]=layer_path;

                cu::xtal::write_poscar(best_report.approximate_tiling_unit_structure,output_dir/tile_path);
                cu::xtal::write_poscar(best_report.approximate_moire_structure,output_dir/root/(lat_to_name(lat)+"_layer.vasp"));

    return;
}

void run_subcommand_twist(const mush::fs::path& input_path, const mush::fs::path& output_dir, const std::vector<double>& angles, int max_lattice_sites, double error_tol, std::string zone, std::string supercells, std::ostream& log)
{
    //GiVe ArGuMenTs LieK aN eDgY tEEn
    std::transform(zone.begin(),zone.end(),zone.begin(),::tolower);
    std::transform(supercells.begin(),supercells.end(),supercells.begin(),::tolower);
    
    
    using LATTICE=mush::MoireStructureReport::LATTICE;
    using ZONE=mush::MoireStructureReport::ZONE;

    ZONE bz= zone=="aligned" ? ZONE::ALIGNED : ZONE::ROTATED;
    if(bz==ZONE::ROTATED)
    {
        assert(zone=="rotated");
    }
                
    /* mush::cautious_create_directory(output_dir); */
    mush::fs::create_directory(output_dir);

    log << "Reading slab from " << input_path << "...\n";
    auto slab = cu::xtal::Structure::from_poscar(input_path);
    //For consistent printing. Should not be necessary for Appriximator classes, which
    //do it internally as well
    mush::make_aligned(&slab);

    mush::json record;
    record["angles"]=angles;
    record["max_lattice_sites"]=max_lattice_sites;
    record["error_tolerance"]=error_tol;


    mush::json twist_record;
    for(const double twist : angles)
    {
        log << "Twising by " << std::fixed << std::setprecision(6) << twist << " degrees ("<<zone<<" Brillouin zone)...\n";
        mush::MoireStructureApproximator moirenator(slab,twist);
        if(max_lattice_sites<moirenator.minimum_lattice_sites(bz))
        {
            log << "Allow minimum possible number of lattice sites in bilayer ("<<moirenator.minimum_lattice_sites(bz)<<")...\n";
        }
        else
        {
            log << "Allow up to "<< max_lattice_sites<<" lattice sites in bilayer...\n";
        }
        moirenator.expand(max_lattice_sites);

        for(auto lat : {LATTICE::ALIGNED,LATTICE::ROTATED})
        {
            if(supercells=="best")
            {

                mush::fs::path root=mush::fs::path(make_twist_dirname(twist));

                mush::fs::create_directories(output_dir/root);
                const auto best_report=moirenator.best_smallest(bz,lat,error_tol);

                auto id=make_twist_id(twist,best_report);
                
                /* twist_record[id][lat_to_name(lat)]=serialize(best_report); */
                /* auto tile_path=root/(lat_to_name(lat)+"_tile.vasp"); */
                /* twist_record[id][lat_to_name(lat)]["tile"]=tile_path; */
                /* auto layer_path=root/(lat_to_name(lat)+"_layer.vasp"); */
                /* twist_record[id][lat_to_name(lat)]["layer"]=layer_path; */

                /* cu::xtal::write_poscar(best_report.approximate_tiling_unit_structure,output_dir/tile_path); */
                /* cu::xtal::write_poscar(best_report.approximate_moire_structure,output_dir/root/(lat_to_name(lat)+"_layer.vasp")); */

                commit_twisted_id(twist,best_report,&twist_record,output_dir,root,lat);

            }

            if(supercells=="all")
            {
                auto best_reports=moirenator.best_of_each_size(bz,lat);
                for(int i=1; i<=best_reports.size(); ++i)
                {
                    auto root=make_target_structure_dir(twist, lat, i);

                    mush::fs::create_directories(output_dir/root);
                    const auto& best_report=best_reports[i-1];

                    /* auto id=make_twist_id(twist,best_report); */

                    /* twist_record[id]=serialize(best_report); */
                    /* auto tile_path=root/(lat_to_name(lat)+"_tile.vasp"); */
                    /* twist_record[id]["tile"]=tile_path; */
                    /* auto layer_path=root/(lat_to_name(lat)+"_layer.vasp"); */
                    /* twist_record[id]["layer"]=layer_path; */

                    /* cu::xtal::write_poscar(best_report.approximate_tiling_unit_structure,output_dir/tile_path); */
                    /* cu::xtal::write_poscar(best_report.approximate_moire_structure,output_dir/root/(lat_to_name(lat)+"_layer.vasp")); */
                    commit_twisted_id(twist,best_report,&twist_record,output_dir,root,lat);
                }
            }
        }
        record["ids"]=twist_record;
    }

    log << "Back up slab structure to " << output_dir / "slab.vasp"
        << "...\n";
    cu::xtal::write_poscar(slab, output_dir / "slab.vasp");

    log << "Save record to "<<output_dir/"record.json"<<"...\n";
    mush::write_json(record,output_dir/"record.json");

    return;
}
