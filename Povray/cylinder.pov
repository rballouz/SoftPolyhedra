#declare tagsam = union{
cylinder {
    <0, 0, -1>,     // Center of one end
    <0, 0, 1>,     // Center of other end
    6            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  
cylinder {
    <0, 0, -1>,     // Center of one end
    <0, 0, 1>,     // Center of other end
    4            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  
cylinder {
    <0, 0, 1>,     // Center of one end
    <0, 0, 10>,     // Center of other end
    0.2            // Radius
    open           // Remove end caps
    pigment { Gray filter 0.5}
  }
  cylinder {
    <0, 0, 1>,     // Center of one end
    <0, 0, 1.0001>,     // Center of other end
    6            // Radius
    pigment { Gray filter 0.5}
  }
  
  cylinder {
    <0, 0, -1>,     // Center of one end
    <0, 0, -0.99999>,     // Center of other end
    6            // Radius
    pigment { Gray filter 0.5}
    }
}