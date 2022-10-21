
using MGroup.MSolve.Discretization.Entities;
using System.Diagnostics;
using MGroup.Constitutive.Structural.Continuum;
using System.Collections.Generic;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.Constitutive.Structural.Transient;
using MGroup.Solvers.Direct;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Logging;

namespace BioApp

{
    class Program
    {
        // define material properties (yperelasticity)
        static double miNormal = 5;//KPa
        static double kappaNormal = 6.667; //Kpa
        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa
        

        static void Main(string[] args)
        {
           
            //Import Mesh 
            var reader = new ComsolMeshReader("../../../meshDataGeom.mphtxt");

            


            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            Comsol3DStaggeredDynamic(reader);

            stopwatch.Stop();
            Console.WriteLine("CPU time : " + stopwatch.ElapsedMilliseconds + "ms");
        }

        public static Model CreateModelFromComsolFile(ComsolMeshReader reader, double miNormal, double kappaNormal, double miTumor, double kappaTumor, Dictionary<int,double> lambda )
        {

            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in reader.NodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }

            //var material = new ConvectionDiffusionProperties(
            //    capacityCoeff: Capacity,
            //    diffusionCoeff: DiffusionCoeff,
            //    convectionCoeff: ConvectionCoeff,
            //    dependentSourceCoeff: DependentProductionCoeff,
            //    independentSourceCoeff: IndependentProductionCoeff);

            var materialNormal = new NeoHookeanMaterial3d(miNormal,kappaNormal);
            var materialTumor = new NeoHookeanMaterial3d(miTumor, kappaTumor);

            var elasticMaterial = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
            var DynamicMaterial = new TransientAnalysisProperties(density: 1, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial,DynamicMaterial);

            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, domainId == 0 ? materialTumor : materialNormal, DynamicMaterial, lambda[elementConnectivity.Key]);
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);

            }

            var faceXYNodes = new List<INode>();
            var faceXZNodes = new List<INode>();
            var faceYZNodes = new List<INode>();

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(0 - node.Z) < 1E-9) faceXYNodes.Add(node);
                if (Math.Abs(0 - node.Y) < 1E-9) faceXZNodes.Add(node);
                if (Math.Abs(0 - node.X) < 1E-9) faceYZNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            foreach (var node in faceXYNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));                
            }
            foreach (var node in faceXZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            }
            foreach (var node in faceYZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            }

            INode maxDistanceNode = null;
            double currentMaxDistance = 0;
            foreach (INode node in model.NodesDictionary.Values)
            {
                double distance = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                if(distance > currentMaxDistance)
                {
                    currentMaxDistance = distance;
                    maxDistanceNode = node;
                }
            }


            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                maxDistanceNode,
                StructuralDof.TranslationX,
                amount: 0.00001
            ));


            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

            return model;
        }

        static void Comsol3DStaggeredDynamic(ComsolMeshReader reader)
        {
            Dictionary<int, double> lambda = new Dictionary<int, double>(reader.ElementConnectivity.Count());
            foreach (var elem in reader.ElementConnectivity)
            {
                lambda.Add(elem.Key, 1);
            }

            var model = new Model[] {
                CreateModelFromComsolFile(reader,  miNormal,  kappaNormal,  miTumor,  kappaTumor,  lambda ),
                //CreateModelFromComsolFile(reader,  miNormal,  kappaNormal,  miTumor,  kappaTumor,  lambda )
            };

            

            var analyzerStates = new GenericAnalyzerState[model.Length];
            var solverFactory = new SkylineSolver.Factory();
            //var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            //var solverFactory = new SkylineSolver.Factory(); //Dense Matrix Solver solves with zero matrices!

            var algebraicModel = new[] {
                solverFactory.BuildAlgebraicModel(model[0]),
                //solverFactory.BuildAlgebraicModel(model[1])
            };

            var solver = new[] {
                solverFactory.BuildSolver(algebraicModel[0]),
                //solverFactory.BuildSolver(algebraicModel[1])
            };

            var problem = new[] {
                new ProblemStructural(model[0], algebraicModel[0], solver[0]),
                //new ProblemConvectionDiffusion(model[1], algebraicModel[1], solver[1])
            };

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(model[0], algebraicModel[0], solver[0], problem[0], numIncrements: 2)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();
            var linearAnalyzers = new[] {
                 new StaticAnalyzer(model[0], algebraicModel[0], solver[0], problem[0], loadControlAnalyzer)
                //new LinearAnalyzer(algebraicModel[1], solver[1], problem[1])
            };


            /*            var analyzer = new[] {
                                        (new NewmarkDynamicAnalyzer.Builder(model[0], algebraicModel[0], solver[0], problem[0], linearAnalyzers[0], timeStep: timeStep, totalTime: totalTime)).Build(),
                                        (new NewmarkDynamicAnalyzer.Builder(model[1], algebraicModel[1], solver[1], problem[1], linearAnalyzers[1], timeStep: timeStep, totalTime: totalTime)).Build(),
                                    };*/

            var analyzer = new[] {
                            (new BDFDynamicAnalyzer.Builder(model[0], algebraicModel[0], solver[0], problem[0], linearAnalyzers[0], timeStep: timeStep, totalTime: totalTime, bdfOrder: 5)).Build(),
                            (new BDFDynamicAnalyzer.Builder(model[1], algebraicModel[1], solver[1], problem[1], linearAnalyzers[1], timeStep: timeStep, totalTime: totalTime, bdfOrder: 5)).Build(),
                        };

            /*            watchDofs = new[] {
                            new List<(INode node, IDofType dof)>(){ (model[0].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), },
                            new List<(INode node, IDofType dof)>(){ (model[1].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), }
                        };*/


            //Sparse Mesh
            watchDofs = new[] {
                new List<(INode node, IDofType dof)>(){ (model[0].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), },
                new List<(INode node, IDofType dof)>(){ (model[1].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), }
            };

            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[1], StructuralDof.TranslationX),
                    (model.NodesDictionary[2], StructuralDof.TranslationZ),

                }, algebraicModel
            );

            linearAnalyzers[0].LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel[0]);
            linearAnalyzers[1].LogFactory = new LinearAnalyzerLogFactory(watchDofs[1], algebraicModel[1]);

            var u1s = new double[(int)(totalTime / timeStep)];
            var u2s = new double[(int)(totalTime / timeStep)];
            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(analyzer, solver, CreateNewModelDynamic, maxStaggeredSteps: 3, tolerance: 1e-5);
            staggeredAnalyzer.Initialize();
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                staggeredAnalyzer.SolveCurrentStep();
                u1s[currentTimeStep] = ((DOFSLog)analyzer[0].ChildAnalyzer.Logs[0]).DOFValues.FirstOrDefault().val;
                u2s[currentTimeStep] = ((DOFSLog)analyzer[1].ChildAnalyzer.Logs[0]).DOFValues.FirstOrDefault().val;

                for (int j = 0; j < analyzer.Length; j++)
                {
                    analyzer[j].AdvanceStep();
                }

                analyzerStates[0] = (analyzer[0] as IParentAnalyzer).CreateState();
                analyzerStates[1] = (analyzer[1] as IParentAnalyzer).CreateState();
            }
        }

        private static void CreateNewModelDynamic(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            double u1 = ((DOFSLog)analyzers[0].ChildAnalyzer.Logs[0]).DOFValues.FirstOrDefault().val;
            double u2 = ((DOFSLog)analyzers[1].ChildAnalyzer.Logs[0]).DOFValues.FirstOrDefault().val;

            /*            model = new[] {
                            Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/OfficialMeshCyprusSparse.mphtxt", ConvectionCoeff: new double[] {1d, 1d, 1d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u2, Capacity: 1d),
                            Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/OfficialMeshCyprusSparse.mphtxt", ConvectionCoeff: new double[] {1d, 1d, 1d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u1, Capacity: 1d)
                        };*/
            //model = new[] {
            //    Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/OfficialMeshCyprusSparse.mphtxt", ConvectionCoeff: new double[] {0d, 0d, 0d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u2, Capacity: 1d),
            //    Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/OfficialMeshCyprusSparse.mphtxt", ConvectionCoeff: new double[] {0d, 0d, 0d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u1, Capacity: 1d)
            //};

            model = new[] {
                Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/3d8Prism.mphtxt", ConvectionCoeff: new double[] {1d, 1d, 1d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u2, Capacity: 1d),
                Comsol3DStaggeredStSt.CreateModelFromComsolFile("../../../Meshes/3d8Prism.mphtxt", ConvectionCoeff: new double[] {1d, 1d, 1d}, DiffusionCoeff: 1d, DependentProductionCoeff: 0d, IndependentProductionCoeff: 1d + u1, Capacity: 1d)
            };

            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            //var solverFactory = new SkylineSolver.Factory();

            var algebraicModel = new[] {
                solverFactory.BuildAlgebraicModel(model[0]),
                solverFactory.BuildAlgebraicModel(model[1])
            };

            solvers[0] = solverFactory.BuildSolver(algebraicModel[0]);
            solvers[1] = solverFactory.BuildSolver(algebraicModel[1]);

            var problem = new[] {
                new ProblemConvectionDiffusion(model[0], algebraicModel[0], solvers[0]),
                new ProblemConvectionDiffusion(model[1], algebraicModel[1], solvers[1])
            };

            var linearAnalyzer = new[] {
                new LinearAnalyzer(algebraicModel[0], solvers[0], problem[0]),
                new LinearAnalyzer(algebraicModel[1], solvers[1], problem[1])
            };

            var oldAnalyzers = analyzers.ToArray();

            //analyzers[0] = (new NewmarkDynamicAnalyzer.Builder(model[0], algebraicModel[0], solvers[0], problem[0], linearAnalyzer[0], timeStep: timeStep, totalTime: totalTime, currentStep: currentTimeStep)).Build();
            //analyzers[1] = (new NewmarkDynamicAnalyzer.Builder(model[1], algebraicModel[1], solvers[1], problem[1], linearAnalyzer[1], timeStep: timeStep, totalTime: totalTime, currentStep: currentTimeStep)).Build();

            analyzers[0] = (new BDFDynamicAnalyzer.Builder(model[0], algebraicModel[0], solvers[0], problem[0], linearAnalyzer[0], timeStep: timeStep, totalTime: totalTime, bdfOrder: 5, currentTimeStep: currentTimeStep)).Build();
            analyzers[1] = (new BDFDynamicAnalyzer.Builder(model[1], algebraicModel[1], solvers[1], problem[1], linearAnalyzer[1], timeStep: timeStep, totalTime: totalTime, bdfOrder: 5, currentTimeStep: currentTimeStep)).Build();

            /*            watchDofs = new[] {
                            new List<(INode node, IDofType dof)>(){ (model[0].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), },
                            n*//*ew List<(INode node, IDofType dof)>(){ (model[1].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), }
                        };*/

            //Sparse Mesh
            watchDofs = new[] {
                new List<(INode node, IDofType dof)>(){ (model[0].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), },
                new List<(INode node, IDofType dof)>(){ (model[1].NodesDictionary[13], ConvectionDiffusionDof.UnknownVariable), }
            };

            linearAnalyzer[0].LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel[0]);
            linearAnalyzer[1].LogFactory = new LinearAnalyzerLogFactory(watchDofs[1], algebraicModel[1]);

            analyzers[0].Initialize(true);
            if (analyzerStates[0] != null)
            {
                analyzers[0].CurrentState = analyzerStates[0];
            }

            analyzers[1].Initialize(true);
            if (analyzerStates[1] != null)
            {
                analyzers[1].CurrentState = analyzerStates[1];
            }
        }
    }
}
