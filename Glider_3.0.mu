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
			"uuid": "E3A7D563-9C87-4E89-AA13-A78F2C39D5A1",
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
			"uuid": "9717D7BA-C7F1-4BC1-8CA7-D22B1517FA15",
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
			"uuid": "CA350524-5D14-4B49-A489-CD32D0BAB14B",
			"type": "WingGeometry",
			"is_main": true,
			"side": "both",
			"span": 7.25,
			"sweep": 0,
			"dihedral": 2.5,
			"mount": 5.3978,
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
			"uuid": "36CAAF8E-39A3-4E16-B244-5A180943023F",
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
				"y": -7.243099606468469,
				"z": -0.316240558398686
			},
			"right_start": {
				"x": 0,
				"y": 7.243099606468469,
				"z": -0.316240558398686
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
			"uuid": "5E2D547B-A20F-49CA-AAE3-942B932DC25F",
			"type": "WingGeometry",
			"is_main": false,
			"side": "both",
			"span": 1,
			"sweep": 0,
			"dihedral": 0,
			"mount": -5.82,
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
			"uuid": "BAB420D8-06AE-4659-B3ED-6404305B7467",
			"type": "WingGeometry",
			"is_main": false,
			"side": "right",
			"span": 0.25,
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
			"uuid": "8818B746-288A-4337-AACB-4BCA0339B45A",
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
			"uuid": "388CFA37-33DB-4950-BAC6-AF18F8F56731",
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
				"uuid": "C752F2B3-F002-41D9-92E6-287F57293ED9",
				"type": "Mesh",
				"name": "Center of Gravity",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-0.10260000079870224,0,-0.014917000196874142,1],
				"geometry": "E3A7D563-9C87-4E89-AA13-A78F2C39D5A1",
				"material": "8818B746-288A-4337-AACB-4BCA0339B45A"
			},
			{
				"uuid": "91F3877D-D508-4D96-8D22-AF7BA567FB5D",
				"type": "Mesh",
				"name": "Aerodynamic Center",
				"visible": false,
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
				"geometry": "9717D7BA-C7F1-4BC1-8CA7-D22B1517FA15",
				"material": "388CFA37-33DB-4950-BAC6-AF18F8F56731"
			},
			{
				"uuid": "C369039A-7997-4B4F-910B-7F3B749E07F2",
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
				"geometry": "CA350524-5D14-4B49-A489-CD32D0BAB14B",
				"material": "9DEEC3B6-DB18-4F96-8E7F-40C207844BAE",
				"children": [
					{
						"uuid": "06237D2F-EB2D-4AD8-86DF-3BADA1778A9C",
						"type": "Mesh",
						"name": "Winglets",
						"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
						"geometry": "36CAAF8E-39A3-4E16-B244-5A180943023F",
						"material": "4C1E6B8A-C503-4F41-8C4E-2E0619E9F625"
					}]
			},
			{
				"uuid": "185BBD51-86BA-44CD-A051-F11F378CEBE5",
				"type": "Mesh",
				"name": "HSTAB",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-3.0569000244140625,0,0,1],
				"geometry": "5E2D547B-A20F-49CA-AAE3-942B932DC25F",
				"material": "9A1A8A8E-B9C4-450D-AF3C-6980AE92C3DA"
			},
			{
				"uuid": "26F7C14D-C6B1-49CD-9E37-41EEF49752EB",
				"type": "Mesh",
				"name": "VSTAB",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-3.394528388977051,0,0,1],
				"geometry": "BAB420D8-06AE-4659-B3ED-6404305B7467",
				"material": "80F8B82B-FC71-4221-A457-8569FA5D54D1"
			}],
		"background": 11184810
	}
}