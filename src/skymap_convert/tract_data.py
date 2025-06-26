"""TODO - Add docstring for tract_data module."""

from dataclasses import dataclass

from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d


@dataclass
class TractData:
    """Represents tract information with basic validation."""

    tract_id: int
    ring: int
    dec_bounds: tuple[float, float]
    ra_bounds: tuple[float, float]
    quad: list[list[float]]
    is_pole: bool | None = None

    def __post_init__(self):
        """Validate tract data after initialization."""
        # Basic validation - declination bounds
        if not (-90 <= self.dec_bounds[0] <= self.dec_bounds[1] <= 90):
            raise ValueError(f"Invalid dec_bounds: {self.dec_bounds}")

        # Validate that we have exactly 4 vertices in quad
        if len(self.quad) != 4:
            raise ValueError(f"Quad must have exactly 4 vertices, got {len(self.quad)}")

        # Validate that each vertex has exactly 2 coordinates
        for i, vertex in enumerate(self.quad):
            if len(vertex) != 2:
                raise ValueError(f"Vertex {i} must have exactly 2 coordinates, got {len(vertex)}")

    def to_convex_polygon(self) -> ConvexPolygon:
        """Convert the tract quad to a ConvexPolygon."""
        unit_vecs = []
        for ra, dec in self.quad:
            unit_vecs.append(UnitVector3d(LonLat.fromDegrees(ra % 360.0, dec)))

        return ConvexPolygon(unit_vecs)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "tract_id": self.tract_id,
            "ring": self.ring,
            "dec_bounds": list(self.dec_bounds),
            "ra_bounds": list(self.ra_bounds),
            "quad": self.quad,
            "is_pole": self.is_pole,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "TractData":
        """Create TractData from dictionary."""
        return cls(
            tract_id=data["tract_id"],
            ring=data["ring"],
            dec_bounds=tuple(data["dec_bounds"]),
            ra_bounds=tuple(data["ra_bounds"]),
            quad=data["quad"],
            is_pole=data.get("is_pole"),
        )
