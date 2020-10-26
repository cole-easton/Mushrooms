﻿using System;
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
	private void buildCap(Func<float, float> curve)
	{
		float curve_0 = curve(0);
		//put the tippy top of the fungus at the origin; we'll 'grow' it downwards 
		vertices.Add(new Vector3(0, 0, 0));
		float subdivisionAngle = Mathf.PI * 2f / axisSubdivisions; //angle between axis subdivisions
		float subdivisionWidth = 1f / capSubdivisions; //width of cap subdivisions
		vertices.Add(subdivisionWidth * new Vector3(1, curve(subdivisionWidth)-curve_0, 0)); //we do the first one out here so that we can make the tris in the for loop too
		for (int i = 1; i < axisSubdivisions; i++)
		{
			float angle = subdivisionAngle * i;
			vertices.Add(subdivisionWidth * new Vector3(Mathf.Cos(angle), curve(subdivisionWidth)-curve_0, Mathf.Sin(angle)));
			tris.AddRange(new int[] { i, 0, i + 1 }); //Unity uses a counter-clockwise winding order, but trig functions also progress counter-clockwise
		}
		tris.AddRange(new int[] { axisSubdivisions, 0, 1 }); // final triangle last to first spoke
		for (int i = 1; i < capSubdivisions; i++)
		{
			int startVertex = vertices.Count;
			vertices.Add((i + 1) * subdivisionWidth * new Vector3(1, curve((i+1)*subdivisionWidth)-curve_0, 0)); //we do the first one out here so that we can make the tris in the for loop too
			for (int j = 1; j < axisSubdivisions; j++)
			{
				float angle = subdivisionAngle * j;
				vertices.Add((i + 1) * subdivisionWidth * new Vector3(Mathf.Cos(angle), curve((i + 1) * subdivisionWidth)-curve_0, Mathf.Sin(angle)));
				int currentVert = 1 + axisSubdivisions * i + j;
				tris.AddRange(new int[] { currentVert - 1, currentVert - axisSubdivisions - 1, currentVert,
					currentVert, currentVert - axisSubdivisions - 1, currentVert - axisSubdivisions});
			}
			tris.AddRange(new int[] { startVertex + axisSubdivisions - 1, startVertex - 1, startVertex,
				startVertex, startVertex-1, startVertex - axisSubdivisions});
			meshFilter.mesh.vertices = vertices.ToArray();
			meshFilter.mesh.SetTriangles(tris.ToArray(), 0);
		}
	}
}