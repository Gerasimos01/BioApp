
using MGroup.MSolve.Discretization.Entities;
using System.Diagnostics;

namespace BioApp

{
    class Program
    {

        static void Main(string[] args)
        {
            //Import Mesh 
            var reader = new ComsolMeshReader("../../../meshDataGeom.mphtxt");

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
           


            stopwatch.Stop();
            Console.WriteLine("CPU time : " + stopwatch.ElapsedMilliseconds + "ms");
        }

        public static Model CreateModelFromComsolFile(ComsolMeshReader reader)
        {
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in reader.NodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }

            var material = new ConvectionDiffusionProperties(
                capacityCoeff: Capacity,
                diffusionCoeff: DiffusionCoeff,
                convectionCoeff: ConvectionCoeff,
                dependentSourceCoeff: DependentProductionCoeff,
                independentSourceCoeff: IndependentProductionCoeff);

            var elementFactory = new ConvectionDiffusionElement3DFactory(material);

            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var element = elementFactory.CreateElement(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2);
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            var topNodes = new List<INode>();
            var bottomNodes = new List<INode>();

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(2 - node.Y) < 1E-9) topNodes.Add(node);
                if (Math.Abs(0 - node.Y) < 1E-9) bottomNodes.Add(node);
            }

            int i = 0;
            var dirichletBCs = new NodalUnknownVariable[topNodes.Count + bottomNodes.Count];
            foreach (var node in topNodes)
            {
                dirichletBCs[i] = new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, 100d);
                i++;
            }
            foreach (var node in bottomNodes)
            {
                dirichletBCs[i] = new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, 50d);
                i++;
            }

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

            return model;
        }

    }
}
