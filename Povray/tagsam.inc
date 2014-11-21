#declare tagsam = union{
cylinder {
    <0, 0, -0.5>,     // Center of one end
    <0, 0, 0.5>,     // Center of other end
    3            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  
cylinder {
    <0, 0, -0.5>,     // Center of one end
    <0, 0, 0.5>,     // Center of other end
    2            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  
cylinder {
    <0, 0, 0.5>,     // Center of one end
    <0, 0, 10>,     // Center of other end
    0.1            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  cylinder {
    <0, 0, 0.5>,     // Center of one end
    <0, 0, 0.5001>,     // Center of other end
    3            // Radius
    pigment { Gray filter 0.5}
  }
  
  cylinder {
    <0, 0, -0.5>,     // Center of one end
    <0, 0, -0.49999>,     // Center of other end
    3            // Radius
    pigment { Gray filter 0.5}
    }
}