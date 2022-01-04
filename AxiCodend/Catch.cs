

namespace AxiCodend
{
    public class Catch
    {
        //=========================
        // variables
        //=========================

        public int[] BlockedMeshes;
        public int Count;


        //=========================
        // constructors
        //=========================

        public Catch() {}

        public Catch(int[] BlockedMeshes)
        {
            this.BlockedMeshes = BlockedMeshes;
            Count = BlockedMeshes.Length;
        }
    }
}
