{
	"MachUp": {
		"version": 4
	},
	"metadata": {
		"version": 4.4,
		"type": "Object",
		"generator": "Object3D.toJSON"
	},
	"geometries": [
		{
			"uuid": "8EF0B389-17DE-419B-A923-69A852863116",
			"type": "SphereBufferGeometry",
			"radius": 0.1,
			"widthSegments": 32,
			"heightSegments": 16,
			"phiStart": 0,
			"phiLength": 6.283185307179586,
			"thetaStart": 0,
			"thetaLength": 3.141592653589793
		},
		{
			"uuid": "D4C4CED4-BADE-4DD0-B2D3-B890F4914961",
			"type": "SphereBufferGeometry",
			"radius": 0.1,
			"widthSegments": 32,
			"heightSegments": 16,
			"phiStart": 0,
			"phiLength": 6.283185307179586,
			"thetaStart": 0,
			"thetaLength": 3.141592653589793
		},
		{
			"uuid": "DCE89460-4E5C-4B28-BB12-71EF3BB8C080",
			"type": "WingGeometry",
			"is_main": true,
			"side": "both",
			"span": 7.25,
			"sweep": 0,
			"dihedral": 0.9049502442,
			"mount": 4.6146,
			"washout": 2.5,
			"root_chord": 0.5,
			"tip_chord": 0.5,
			"root_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.0455,
						"CL_alpha": 5.736,
						"Cm_L0": -0.0388,
						"Cm_alpha": -0.08,
						"CD0": 0.0064,
						"CD0_L": -0.0021,
						"CD0_L2": 0.0062,
						"CL_max": 1.4
					}
				}
			},
			"tip_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.0455,
						"CL_alpha": 5.736,
						"Cm_L0": -0.0388,
						"Cm_alpha": -0.08,
						"CD0": 0.0064,
						"CD0_L": -0.0021,
						"CD0_L2": 0.0062,
						"CL_max": 1.4
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"right_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		},
		{
			"uuid": "6E862146-4050-428A-B110-891432399A50",
			"type": "WingGeometry",
			"is_main": true,
			"side": "both",
			"span": 0.5,
			"sweep": 0,
			"dihedral": 90,
			"mount": 0,
			"washout": 0,
			"root_chord": 0.5,
			"tip_chord": 0.5,
			"root_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.0455,
						"CL_alpha": 5.736,
						"Cm_L0": -0.0388,
						"Cm_alpha": -0.08,
						"CD0": 0.0064,
						"CD0_L": -0.0021,
						"CD0_L2": 0.0062,
						"CL_max": 1.4
					}
				}
			},
			"tip_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.0455,
						"CL_alpha": 5.736,
						"Cm_L0": -0.0388,
						"Cm_alpha": -0.08,
						"CD0": 0.0064,
						"CD0_L": -0.0021,
						"CD0_L2": 0.0062,
						"CL_max": 1.4
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": -7.249095719593477,
				"z": -0.11450435874470508
			},
			"right_start": {
				"x": 0,
				"y": 7.249095719593477,
				"z": -0.11450435874470508
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		},
		{
			"uuid": "CD3242B1-4DBA-46BE-B99A-15537A09DCB9",
			"type": "WingGeometry",
			"is_main": false,
			"side": "both",
			"span": 1,
			"sweep": 0,
			"dihedral": 0,
			"mount": 0,
			"washout": 0,
			"root_chord": 0.5,
			"tip_chord": 0.35,
			"root_airfoil": {
				"NACA 0012": {
					"properties": {
						"type": "linear",
						"alpha_L0": 0,
						"CL_alpha": 6.118,
						"Cm_L0": 0,
						"Cm_alpha": 0,
						"CD0": 0.0058,
						"CD0_L": 0,
						"CD0_L2": 0.0059,
						"CL_max": 1.2
					}
				}
			},
			"tip_airfoil": {
				"NACA 0012": {
					"properties": {
						"type": "linear",
						"alpha_L0": 0,
						"CL_alpha": 6.118,
						"Cm_L0": 0,
						"Cm_alpha": 0,
						"CD0": 0.0058,
						"CD0_L": 0,
						"CD0_L2": 0.0059,
						"CL_max": 1.2
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"right_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		},
		{
			"uuid": "2739DF8E-B9A3-4A97-88EC-07D19AABAE3D",
			"type": "WingGeometry",
			"is_main": false,
			"side": "right",
			"span": 0.75,
			"sweep": 0,
			"dihedral": 90,
			"mount": 0,
			"washout": 0,
			"root_chord": 0.5,
			"tip_chord": 0.35,
			"root_airfoil": {
				"NACA 0012": {
					"properties": {
						"type": "linear",
						"alpha_L0": 0,
						"CL_alpha": 6.118,
						"Cm_L0": 0,
						"Cm_alpha": 0,
						"CD0": 0.0058,
						"CD0_L": 0,
						"CD0_L2": 0.0059,
						"CL_max": 1.2
					}
				}
			},
			"tip_airfoil": {
				"NACA 0012": {
					"properties": {
						"type": "linear",
						"alpha_L0": 0,
						"CL_alpha": 6.118,
						"Cm_L0": 0,
						"Cm_alpha": 0,
						"CD0": 0.0058,
						"CD0_L": 0,
						"CD0_L2": 0.0059,
						"CL_max": 1.2
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"right_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		}],
	"materials": [
		{
			"uuid": "53628DD4-DC3B-40AC-ABFF-10499FFFB764",
			"type": "MeshStandardMaterial",
			"color": 16711680,
			"roughness": 0.5,
			"metalness": 0.5,
			"emissive": 16711680,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "652B0758-568C-420E-B260-685E20C25E00",
			"type": "MeshStandardMaterial",
			"color": 6684927,
			"roughness": 0.5,
			"metalness": 0.5,
			"emissive": 6684927,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "9DEEC3B6-DB18-4F96-8E7F-40C207844BAE",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "4C1E6B8A-C503-4F41-8C4E-2E0619E9F625",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "9A1A8A8E-B9C4-450D-AF3C-6980AE92C3DA",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "80F8B82B-FC71-4221-A457-8569FA5D54D1",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		}],
	"object": {
		"uuid": "BA245872-5044-46C3-A4C6-59C7D2F11880",
		"type": "Origin",
		"name": "MyAirplane",
		"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
		"children": [
			{
				"uuid": "C554EA28-85B0-4A97-9555-F25A70EBC2DB",
				"type": "Mesh",
				"name": "Center of Gravity",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0.10000000149011612,0,-0.015443000011146069,1],
				"geometry": "8EF0B389-17DE-419B-A923-69A852863116",
				"material": "53628DD4-DC3B-40AC-ABFF-10499FFFB764"
			},
			{
				"uuid": "BBC85F9B-94E1-433C-A11E-B7D11128FBDB",
				"type": "Mesh",
				"name": "Aerodynamic Center",
				"visible": false,
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
				"geometry": "D4C4CED4-BADE-4DD0-B2D3-B890F4914961",
				"material": "652B0758-568C-420E-B260-685E20C25E00"
			},
			{
				"uuid": "4FCF1C3F-3BC0-40E9-B818-FF5BCF1920FE",
				"type": "PointLight",
				"name": "PointLight",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,10,10,-20,1],
				"color": 16777215,
				"intensity": 1,
				"distance": 0,
				"decay": 1
			},
			{
				"uuid": "4B0FE6F9-7632-49F9-8961-828A2CC7F645",
				"type": "Mesh",
				"name": "Main_Wing",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
				"geometry": "DCE89460-4E5C-4B28-BB12-71EF3BB8C080",
				"material": "9DEEC3B6-DB18-4F96-8E7F-40C207844BAE",
				"children": [
					{
						"uuid": "06237D2F-EB2D-4AD8-86DF-3BADA1778A9C",
						"type": "Mesh",
						"name": "Winglets",
						"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
						"geometry": "6E862146-4050-428A-B110-891432399A50",
						"material": "4C1E6B8A-C503-4F41-8C4E-2E0619E9F625"
					}]
			},
			{
				"uuid": "185BBD51-86BA-44CD-A051-F11F378CEBE5",
				"type": "Mesh",
				"name": "HSTAB",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-3,0,0,1],
				"geometry": "CD3242B1-4DBA-46BE-B99A-15537A09DCB9",
				"material": "9A1A8A8E-B9C4-450D-AF3C-6980AE92C3DA"
			},
			{
				"uuid": "26F7C14D-C6B1-49CD-9E37-41EEF49752EB",
				"type": "Mesh",
				"name": "VSTAB",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-3.394528388977051,0,0,1],
				"geometry": "2739DF8E-B9A3-4A97-88EC-07D19AABAE3D",
				"material": "80F8B82B-FC71-4221-A457-8569FA5D54D1"
			}],
		"background": 11184810
	}
}