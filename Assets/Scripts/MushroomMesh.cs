using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MushroomMesh : MonoBehaviour
{
	public int capSubdivisions = 4;
	public int axisSubdivisions = 12;
	
	private List<Vector3> vertices = new List<Vector3>();
	private List<int> tris = new List<int>();
	private MeshFilter meshFilter;

    // Start is called before the first frame update
    void Start()
    {
		meshFilter = gameObject.GetComponent<MeshFilter>();
		buildCap(x => -x*x);
		
		meshFilter.mesh.RecalculateNormals();  
    }

    // Update is called once per frame
    void Update()
    {
        
    }

	/// <summary>
	/// Builds the top of the mushroom's cap and adds the vertices and trianges to the submesh at index 0.  Does not recalulate normals
	/// precondition: vertices.Count == 0
	/// </summary>
	/// <param name="curve">The function, evaluated over [0, 1] that defines the shape of the cap.</param>
	/// <param name = "includeGradient">Whether the edge of the cap should round around to the bottom</param>
	private void buildCap(Func<float, float> curve, bool includeGradient = true)
	{
		float curve_0 = curve(0);
		//put the tippy top of the fungus at the origin; we'll 'grow' it downwards 
		vertices.Add(new Vector3(0, 0, 0));
		float subdivisionAngle = Mathf.PI * 2f / axisSubdivisions; //angle between axis subdivisions
		float subdivisionWidth = 1f / capSubdivisions; //width of cap subdivisions
		vertices.Add(new Vector3(subdivisionWidth, curve(subdivisionWidth)-curve_0, 0)); //we do the first one out here so that we can make the tris in the for loop too
		for (int i = 1; i < axisSubdivisions; i++)
		{
			float angle = subdivisionAngle * i;
			vertices.Add(new Vector3(subdivisionWidth * Mathf.Cos(angle), curve(subdivisionWidth)-curve_0, subdivisionWidth * Mathf.Sin(angle)));
			tris.AddRange(new int[] { i, 0, i + 1 }); //Unity uses a counter-clockwise winding order, but trig functions also progress counter-clockwise
		}
		tris.AddRange(new int[] { axisSubdivisions, 0, 1 }); // final triangle last to first spoke
		for (int i = 1; i < capSubdivisions; i++)
		{
			int startVertex = vertices.Count;
			vertices.Add(new Vector3((i + 1) * subdivisionWidth, curve((i+1)*subdivisionWidth)-curve_0, 0)); //we do the first one out here so that we can make the tris in the for loop too
			for (int j = 1; j < axisSubdivisions; j++)
			{
				float angle = subdivisionAngle * j;
				vertices.Add(new Vector3((i + 1) * subdivisionWidth * Mathf.Cos(angle), curve((i + 1) * subdivisionWidth)-curve_0, (i + 1) * subdivisionWidth * Mathf.Sin(angle)));
				int currentVert = 1 + axisSubdivisions * i + j;
				tris.AddRange(new int[] { currentVert - 1, currentVert - axisSubdivisions - 1, currentVert,
					currentVert, currentVert - axisSubdivisions - 1, currentVert - axisSubdivisions});
			}
			tris.AddRange(new int[] { startVertex + axisSubdivisions - 1, startVertex - 1, startVertex,
				startVertex, startVertex-1, startVertex - axisSubdivisions});
		}
		if (includeGradient)
		{
			float direction = (Mathf.Atan2(curve(1.001f) - curve(0.999f), 0.002f) + Mathf.PI*2)%(Mathf.PI*2); //approximation of angle of tangent 
			Debug.Log(direction);
			float radius = 1, y = curve(1);
			int i = capSubdivisions;
			while (direction > Mathf.PI / 2)
			{
				radius += Mathf.Cos(direction) / 20;
				y += Mathf.Sin(direction) / 20;
				int startVertex = vertices.Count;
				vertices.Add(new Vector3(radius, y, 0)); //we do the first one out here so that we can make the tris in the for loop too
				for (int j = 1; j < axisSubdivisions; j++)
				{
					float angle = subdivisionAngle * j;
					vertices.Add(new Vector3(radius * Mathf.Cos(angle), y, radius*Mathf.Sin(angle)));
					int currentVert = 1 + axisSubdivisions * i + j;
					tris.AddRange(new int[] { currentVert - 1, currentVert - axisSubdivisions - 1, currentVert,
					currentVert, currentVert - axisSubdivisions - 1, currentVert - axisSubdivisions});
				}
				tris.AddRange(new int[] { startVertex + axisSubdivisions - 1, startVertex - 1, startVertex,
				startVertex, startVertex-1, startVertex - axisSubdivisions});
				direction -= Mathf.PI / 10;
				i++;
			}

		}
		meshFilter.mesh.vertices = vertices.ToArray();
		meshFilter.mesh.SetTriangles(tris.ToArray(), 0);
	}

	private void buildGills(Func<float, float> curve, float length)
	{
		int currentVert = vertices.Count;

	}
}
