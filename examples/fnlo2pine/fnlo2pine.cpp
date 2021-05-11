#include <cstdlib>
#include <string>

#include <fastnlotk/fastNLOLHAPDF.h>

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        return EXIT_FAILURE;
    }

    std::string in(argv[1]);
    std::string out(argv[2]);

    auto file = fastNLOLHAPDF(in, "NNPDF31_nlo_as_0118_luxqed");
    //auto ref_table = file.GetReferenceTable(fastNLO::kNextToNextToLeading);
    auto coef_table = file.GetCoeffTable(0);

    assert(coef_table != nullptr);

    auto converted = dynamic_cast <fastNLOCoeffAddFix*> (coef_table);

    assert(converted != nullptr);

    auto descriptions = coef_table->GetContributionDescription();

    for (auto const& desc : descriptions)
    {
        std::cout << "-- " << desc << '\n';
    }
}
