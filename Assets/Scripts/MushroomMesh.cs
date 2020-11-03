using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

public class MushroomMesh : MonoBehaviour
{
	public int capSubdivisions = 4;
	public int axisSubdivisions = 12;
	public int stipeSubdivisionDensity;

	//since Func cannot be shown in the inspector, all fields that control the shape/form of the mushroom are hidden from it
	[HideInInspector]
	public Func<float, float> capCurve = x => -x*x;
	[HideInInspector]
	public bool roundedCapEdge = true;
	[HideInInspector]
	public float gillOffset = -0.2f;
	[HideInInspector]
	public float gillInnerRadius = 0.35f;
	[HideInInspector]
	/// <summary>
	/// the stipe curve, evaluated over [1, length] that defines the shape of the stipe. the x and z output control the position, 
	/// and the y controls the radius (relative to the initial radius).  The function will be automatically translated so that curve(0) == new Vector3(0, 1, 0)
	/// </summary>
	public Func<float, Vector3> stipeCurve = y => new Vector3(0.2f * Mathf.Cos(1.5f * y), -0.5f * (y - 1) * (y - 1) + 1.5f, 0.2f * Mathf.Sin(1.5f * y));
	[HideInInspector]
	public float stipeLength = 1.5f;


	private List<Vector3> vertices = new List<Vector3>();
	private List<Vector2> uVs = new List<Vector2>();
	private MeshFilter meshFilter;
	private MeshRenderer meshRenderer;

    // Start is called before the first frame update
    void Start()
    {
		RandomizeParams();
		meshFilter = gameObject.GetComponent<MeshFilter>();
		meshFilter.mesh.subMeshCount = 3;

		meshRenderer = gameObject.GetComponent<MeshRenderer>();

		BuildCap(capCurve, roundedCapEdge);
		BuildGills(gillInnerRadius, gillOffset);
		BuildStipe(stipeCurve, stipeLength);
		CorrectPosition();
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
	private void BuildCap(Func<float, float> curve, bool includeGradient = true)
	{
		List<int> tris = new List<int>();
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
			float deltaY = curve(1.001f) - curve(0.999f);
			float direction;
			if (deltaY == 0)
			{
				direction = 1.5f * Mathf.PI;
			}
			else
			{
				direction = (Mathf.Atan2(curve(1.001f) - curve(0.999f), 0.002f) + Mathf.PI * 2) % (Mathf.PI * 2); //approximation of angle of tangent 
			}
			float radius = 1, y = curve(1)-curve_0;
			int i = capSubdivisions;
			while (direction > Mathf.PI / 2)
			{
				radius += Mathf.Cos(direction) / 40; //fineness (influences size too)
				y += Mathf.Sin(direction) / 40;
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
				direction -= Mathf.PI / 10; // sharpness of curvature (sharper means smaller, too)
				i++;
			}

		}
		meshFilter.mesh.vertices = vertices.ToArray();
		meshFilter.mesh.SetTriangles(tris.ToArray(), 0);

		//build UVs
		//calculare arc length
		float arcLength = 0;
		for (int i = axisSubdivisions; i < vertices.Count; i += axisSubdivisions)
		{
			arcLength += Vector3.Distance(vertices[i], vertices[i - axisSubdivisions]);
		}
		uVs.Add(new Vector2(0.5f, 0.5f)); //center point of cap is in the middle
		float radMultiplier = 0.49f / arcLength;
		float currentRad = 0;
		for (int i = 1; i < vertices.Count; i++)
		{
			int lastIndex = i < axisSubdivisions ? 0 : i - axisSubdivisions;
			if ((i-1)%axisSubdivisions == 0)
				currentRad += Vector3.Distance(vertices[i], vertices[lastIndex]);
			Vector2 uv = new Vector2(vertices[i].x, vertices[i].z);
			uVs.Add((uv.normalized * radMultiplier * currentRad) + new Vector2(0.5f, 0.5f));
		}
		meshFilter.mesh.uv = uVs.ToArray();
	}

	private void BuildGills(float radius, float yOffset = 0)
	{
		vertices.AddRange(vertices.GetRange(vertices.Count - axisSubdivisions, axisSubdivisions)); // duplicate ring for hard edge
		float uvRadMultiplier = .5f/Mathf.Sqrt(vertices[vertices.Count - 1].x * vertices[vertices.Count - 1].x + vertices[vertices.Count - 1].z * vertices[vertices.Count - 1].z);
		float subdivisionAngle = Mathf.PI * 2f / axisSubdivisions; //angle between axis subdivisions
		for (int i = 0; i < axisSubdivisions; i++)
		{
			uVs.Add(new Vector2(.5f * Mathf.Cos(subdivisionAngle * i) + 0.5f, .5f * Mathf.Sin(subdivisionAngle * i) + 0.5f));
		}
		List<int> tris = new List<int>();
		int startVert = vertices.Count;
		float y = vertices[vertices.Count - 1].y + yOffset;
		vertices.Add(new Vector3(radius, y, 0));
		uVs.Add(new Vector2(radius * uvRadMultiplier + 0.5f, 0.5f));
		for (int j = 1; j < axisSubdivisions; j++)
		{
			float angle = subdivisionAngle * j;
			vertices.Add(new Vector3(radius * Mathf.Cos(angle), y, radius * Mathf.Sin(angle)));
			uVs.Add(new Vector2(radius*uvRadMultiplier*Mathf.Cos(angle) + 0.5f, radius*uvRadMultiplier*Mathf.Sin(angle) + 0.5f));
			int currentVert = startVert + j;
			tris.AddRange(new int[] { currentVert - 1, currentVert - axisSubdivisions - 1, currentVert,
					currentVert, currentVert - axisSubdivisions - 1, currentVert - axisSubdivisions});
		}
		tris.AddRange(new int[] { startVert + axisSubdivisions - 1, startVert - 1, startVert,
				startVert, startVert-1, startVert - axisSubdivisions});
		meshFilter.mesh.vertices = vertices.ToArray(); 
		meshFilter.mesh.uv = uVs.ToArray();
		meshFilter.mesh.SetTriangles(tris.ToArray(), 1);

	}

	/// <summary>
	/// builds the stipe (stem) of the mushroom and calculates the normals of the entire mesh.  Do not use RecalculateNormals after runnung this method, 
	/// as it will cause undesirable hard edges.
	/// </summary>
	/// <param name="curve">the curve, evaluated over [1, length] that defines the shape of the stipe. the x and z output control the position, 
	/// and the y controls the radius (relative to the initial radius).  The function will be automatically translated so that curve(0) == new Vector3(0, 1, 0) </param>
	/// <param name="length">The height if the stipe (not the arc length)</param>
	private void BuildStipe(Func<float, Vector3> curve, float length = 1) //the stipe is the "stem" of the mushroom
	{
		vertices.AddRange(vertices.GetRange(vertices.Count - axisSubdivisions, axisSubdivisions)); // duplicate ring so we can assign unique UV coords
		int firstSeamIndex = vertices.Count - axisSubdivisions;
		vertices.Add(vertices[firstSeamIndex]);
		for (int j = 0; j <= axisSubdivisions; j++)
		{
			uVs.Add(new Vector2((float)j / axisSubdivisions, 0));
		}
		List<int> tris = new List<int>();
		float subdivisionAngle = Mathf.PI * 2f / axisSubdivisions; //angle between axis subdivisions
		float subdivisionHeight = 1f / stipeSubdivisionDensity;
		float y_0 = vertices[vertices.Count - 1].y;
		float rad_0 = Mathf.Sqrt(vertices[vertices.Count-1].x * vertices[vertices.Count-1].x + vertices[vertices.Count-1].z * vertices[vertices.Count-1].z);
		float radOffset = 1 - curve(0).y;
		float xOffset = -curve(0).x;
		float zOffset = -curve(0).z;
		int subdivisions = (int)(length * stipeSubdivisionDensity + 0.5f);
		for (int i = 1; i < subdivisions; i++)
		{
			Vector3 currentParams = curve(subdivisionHeight * i);
			float radius = currentParams.y;
			int startVert = vertices.Count;
			vertices.Add(new Vector3(rad_0*(radius+radOffset) + currentParams.x + xOffset, y_0 - subdivisionHeight*i, currentParams.z + zOffset));
			uVs.Add(new Vector2(0, (float)i/(subdivisions-1)));
			for (int j = 1; j < axisSubdivisions + 1; j++)
			{
				float angle = subdivisionAngle * j; 
				vertices.Add(new Vector3(rad_0 * (radius + radOffset) * Mathf.Cos(angle) + currentParams.x + xOffset, y_0-subdivisionHeight*i, rad_0 * (radius + radOffset) * Mathf.Sin(angle) + currentParams.z + zOffset));
				uVs.Add(new Vector2((float)j / axisSubdivisions, (float)i / (subdivisions-1)));
				int currentVert = startVert + j;
				
				tris.AddRange(new int[] { currentVert - 1, currentVert - axisSubdivisions - 2, currentVert,
					currentVert, currentVert - axisSubdivisions - 2, currentVert - axisSubdivisions - 1});
			}
		}
		meshFilter.mesh.vertices = vertices.ToArray();
		meshFilter.mesh.uv = uVs.ToArray();
		meshFilter.mesh.SetTriangles(tris.ToArray(), 2);
		meshFilter.mesh.RecalculateNormals();
		//smooth out hard edge
		Vector3[] normals = meshFilter.mesh.normals; // "To make changes to the normals it is important to copy the normals from the Mesh." --https://docs.unity3d.com/ScriptReference/Mesh-normals.html 11/2/2020
		for (int i = firstSeamIndex; i + axisSubdivisions < vertices.Count; i += axisSubdivisions + 1)
		{
			Vector3 avgNorm = (normals[i] + normals[i + axisSubdivisions]) / 2f;
			normals[i] = normals[i + axisSubdivisions] = avgNorm;
		}
		meshFilter.mesh.normals = normals;
	}

	/// <summary>
	/// corrects the position of the mushroom so that that bottom of the stipe touches the ground
	/// </summary>
	private void CorrectPosition()
	{
		float height = capCurve(0) - capCurve(1) - gillOffset + stipeLength;
		transform.Translate(height * Vector3.up);
	}

	/// <summary>
	/// Randomly generates new parameters for the shape/form of the mushroom using values and distrubutions according 
	/// to the author's (Cole R. Easton) subjective aesthetic opinions
	/// </summary>
	public void RandomizeParams()
	{
		float rng = Random.value;
		if (rng < 0.25f) //concave cap
		{
			Func<float, float>[] funcs = new Func<float, float>[] {x=>0.5f*x*x, x=>Mathf.Sqrt(x*x+1), x=>-Mathf.Sqrt(2-x*x),
				x =>Mathf.Pow(x, 4),  x=>0.25f*Mathf.Pow(x, 4), x => 1.92901f*x*x*x - 2.31481f*x*x - 0.123457f };
			capCurve = funcs[Random.Range(0, funcs.Length)];
			roundedCapEdge = false;
			float dydx = (capCurve(1.001f) - capCurve(0.999f)) / 0.002f;
			gillInnerRadius = Random.Range(0.1f, 0.5f);
			float minOffset = -dydx * (1 - gillInnerRadius);
			gillOffset = Random.Range(minOffset, minOffset * 1.2f);
		}
		else //convex cap
		{
			float quadraticCoeff = -Gaussian(.8f, 0.2f);
			float quarticCoeff = -Mathf.Pow(2, Random.Range(-2f, 0f));
			float circleRad = Random.Range(1f, 2.5f);
			float bellSteepness = Random.Range(0.5f, 4f);
			float cosScalar = Random.Range(0.8f, Mathf.PI);
			float secCoeff = Random.Range(0.5f, 1.2f);
			Func<float, float>[] funcs = new Func<float, float>[] {x=>quadraticCoeff*x*x,  x=>quarticCoeff*Mathf.Pow(x, 4), x=>Mathf.Sqrt(circleRad-x*x), 
				  x => Mathf.Exp(-bellSteepness*x*x), x => Mathf.Cos(cosScalar*x), x => -secCoeff/Mathf.Cos(x) };
			capCurve = funcs[Random.Range(0, funcs.Length)];
			roundedCapEdge = Random.value < 0.8f;
			gillInnerRadius = Mathf.Pow(Gaussian(0.59f, 0.12f), 2);
			if (gillInnerRadius < 0.05f)
				gillInnerRadius = 0.05f;
			if (gillInnerRadius > 0.95f)
				gillInnerRadius = 0.95f;
			gillOffset = Random.Range(-0.4f, 0.4f);
			if (gillOffset >= capCurve(gillInnerRadius)-capCurve(1)) //so we don't get stipes sticking throught the cap
			{
				gillOffset = capCurve(gillInnerRadius)-capCurve(1) - 0.01f;
			}
		}
		stipeLength = Gaussian(1.5f, 0.5f);
		float cosOffset = Random.Range(0, 2 * Mathf.PI);
		float sinOfffset = Random.Range(0, 2 * Mathf.PI);
		float val = Random.value;
		if (val < 0.4f) //convex
		{
			float fatSpot = Random.Range(0.5f, 0.75f); // 2/3 in default curve
			stipeCurve = y => new Vector3(0.2f * Mathf.Cos(1.5f * y + cosOffset), -0.5f * Mathf.Pow(y - stipeLength * fatSpot, 2f) + 1.5f, 0.2f * Mathf.Sin(1.5f * y + sinOfffset));
		}
		else if (val < 0.7) //concave
		{
			float verticalOffset = Random.Range(0.25f, 0.5f); //.5 is symmetrical; .25 is quite bottom heavy
			float subtleness = .1f * Mathf.Pow(2f, Random.Range(-2f, 1f));
			stipeCurve = y => new Vector3(0.2f * Mathf.Cos(1.5f * y + cosOffset), Mathf.Pow(y/stipeLength - verticalOffset, 4f) + subtleness, 0.2f * Mathf.Sin(1.5f * y + sinOfffset));
		}
		else 
		{
			stipeCurve = y => new Vector3(0.2f * Mathf.Cos(1.5f * y + cosOffset), 1f, 0.2f * Mathf.Sin(1.5f * y + sinOfffset));
		}
	}

	private float Gaussian(float mean, float stdDev)
	{
		float val1 = Random.Range(0f, 1f);
		float val2 = Random.Range(0f, 1f);
		float gaussValue = Mathf.Sqrt(-2.0f * Mathf.Log(val1)) * Mathf.Sin(2.0f * Mathf.PI * val2);
		return mean + stdDev * gaussValue;
	}
}
