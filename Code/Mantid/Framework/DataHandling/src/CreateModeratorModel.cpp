/*WIKI*
Creates a model of a moderator using the given parameters. The parameters are given as a string to allow flexibility for each moderator model having different parameterisation.

Available models with parameter names:
* IkedaCarpenterModerator -
** TiltAngle - The tilt angle in degrees;
** TauF - The fast decay coefficient in microseconds
** TauS - The slow decay coefficient in microseconds
** R - The mixing coefficient
*WIKI*/

#include "MantidDataHandling/CreateModeratorModel.h"

#include "MantidAPI/IkedaCarpenterModerator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/MandatoryValidator.h"

namespace Mantid
{
namespace DataHandling
{
  using API::WorkspaceProperty;
  using API::MatrixWorkspace_sptr;
  using Kernel::Direction;
  using Kernel::ListValidator;
  using Kernel::MandatoryValidator;

  // Register the algorithm into the AlgorithmFactory
  DECLARE_ALGORITHM(CreateModeratorModel)
  
  //----------------------------------------------------------------------------------------------
  /// Algorithm's name for identification. @see Algorithm::name
  const std::string CreateModeratorModel::name() const { return "CreateModeratorModel";};
  
  /// Algorithm's version for identification. @see Algorithm::version
  int CreateModeratorModel::version() const { return 1;};
  
  /// Algorithm's category for identification. @see Algorithm::category
  const std::string CreateModeratorModel::category() const { return "DataHandling";}

  //----------------------------------------------------------------------------------------------
  /// Sets documentation strings for this algorithm
  void CreateModeratorModel::initDocs()
  {
    this->setWikiSummary("Creates the given moderator model and attaches it to the input workspace.");
    this->setOptionalMessage("Creates the given moderator model and attaches it to the input workspace.");
  }

  //----------------------------------------------------------------------------------------------
  /** Initialize the algorithm's properties.
   */
  void CreateModeratorModel::init()
  {
    declareProperty(new WorkspaceProperty<>("Workspace","",Direction::InOut), "An input workspace.");

    std::vector<std::string> keys(1, "IkedaCarpenterModerator");
    declareProperty("ModelType", "", boost::make_shared<ListValidator<std::string> >(keys),
                    "The string identifier for the model", Direction::Input);

    declareProperty("Parameters", "", boost::make_shared<MandatoryValidator<std::string> >(),
                    "The parameters for the model as comma-separated list of name=value pairs");
  }

  //----------------------------------------------------------------------------------------------
  /** Execute the algorithm.
   */
  void CreateModeratorModel::exec()
  {
    const std::string modelType = getProperty("ModelType");
    if(modelType != "IkedaCarpenterModerator")
    {
      throw std::invalid_argument("Invalid moderator model type.");
    }

    API::ModeratorModel *moderator = new API::IkedaCarpenterModerator;
    moderator->initialize(getProperty("Parameters"));

    MatrixWorkspace_sptr workspace = getProperty("Workspace");
    workspace->setModeratorModel(moderator);
  }



} // namespace DataHandling
} // namespace Mantid
